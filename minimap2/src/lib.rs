//! High-level, safe Rust wrapper for minimap2

use anyhow::{bail, Result};
use noodles_fastq as fastq;
use noodles_sam::{
    alignment::RecordBuf, 
    header::Header,
    header::record::value::{map::ReferenceSequence, Map},
    alignment::record_buf::Cigar,
    alignment::record::cigar::{op::Kind, Op},
};
use std::{
    ffi::{c_char, c_int, CString},
    path::{Path, PathBuf},
    sync::Arc,
    cmp::{max, min},
    num::NonZeroUsize,
    collections::BTreeMap,
};

pub use minimap2_sys;

/// Main interface for performing read alignments using minimap2
#[derive(Debug)]
pub struct Minimap2Aligner {
    index: Arc<Minimap2Index>,
    opts: Minimap2Opts,
    inner: InnerAligner,
}

impl Minimap2Aligner {
    /// Loads the genome index into memory.
    ///
    /// # Arguments
    /// - `opts` - Options controlling alignment parameters.
    ///
    /// # Returns
    /// - A `Minimap2Aligner` object initialized with the provided options.
    pub fn new(opts: Minimap2Opts) -> Result<Self> {
        let index = Minimap2Index::load(&opts)?;
        let inner = InnerAligner::new(&opts)?;
        Ok(Self {
            index: Arc::new(index),
            opts,
            inner,
        })
    }
    
    /// Create a new aligner with specified thread count for optimal performance
    /// Currently, extra_params has not been implemented
    pub fn with_threads(opts: Minimap2Opts, num_threads: usize) -> Result<Self> {
        let mut aligner = Self::new(opts)?;

        // Store thread count for batch operations
        // --threads=num
        aligner.opts.extra_params.push(format!("--threads={}", num_threads));
        Ok(aligner)
    }


    /// Retrieves the SAM header information.
    ///
    /// # Returns
    /// - A reference to the SAM header containing reference sequence metadata
    pub fn get_header(&self) -> &Header {
        &self.index.header
    }

    /// Align long read assay (single-end)
    /// Main function for aligning a single read
    pub fn align_read(&mut self, fq: &fastq::Record) -> Result<Vec<RecordBuf>> {
        self.inner.align_read(&self.index, fq)
    }
    
    /// Batch alignment of multiple fastq query records
    /// Every thread gets its own InnerAligner with independent mm_tbuf_t (provided in minimap2_sys)
    /// Chunks are separated units of batch fastq records
    pub fn align_batch(&self, records: &[fastq::Record], num_threads: usize) -> Result<Vec<Vec<RecordBuf>>> {
        if records.is_empty() {
            return Ok(Vec::new());
        }
        // Use num_threads, TODO: consider using Rayon for parallel processing
        let num_threads = if num_threads == 0 {
            let available = std::thread::available_parallelism()?.get();
            max(1, min(available, 16))
        } else {
            max(1, min(num_threads, 32))
        };
        
        // Use scoped threads to process batches in parallel
        let chunk_size = (records.len() + num_threads - 1) / num_threads;
        let chunks: Vec<_> = records.chunks(chunk_size).collect();
        
        let results = std::thread::scope(|s| {
            let handles: Vec<_> = chunks
                .into_iter()
                .map(|chunk| {
                    let index = self.index.clone();
                    let opts = self.opts.clone();
                    // Each thread gets its own InnerAligner with independent mm_tbuf_t
                    s.spawn(move || {
                        let mut aligner = InnerAligner::new(&opts)?;
                        let mut thread_results = Vec::new();
                        for record in chunk {
                            let alignments = aligner.align_read(&index, record)?;
                            thread_results.push(alignments);
                        }
                        Ok::<Vec<Vec<RecordBuf>>, anyhow::Error>(thread_results)
                    })
                })
                .collect();
            
            let mut all_results = Vec::new();
            for handle in handles {
                match handle.join() {
                    Ok(thread_result) => all_results.extend(thread_result?),
                    Err(_) => bail!("Thread panicked during batch alignment"),
                }
            }
            Ok::<Vec<Vec<RecordBuf>>, anyhow::Error>(all_results)
        })?;
        Ok(results)
    }
    
    }

/// This will return a copy of the aligner with the same reference index.
/// Returning a new aligner is necessary for multithreaded alignment, as the
/// underlying foreign aligner is not thread-safe.
impl Clone for Minimap2Aligner {
    fn clone(&self) -> Self {
        Self {
            index: self.index.clone(),
            opts: self.opts.clone(),
            inner: InnerAligner::new(&self.opts)
                .expect("Failed to create inner aligner during clone"),
        }
    }
}

/// minimap2 presets for DNA & mRNA/cDNA
/// 
/// - Oxford Nanopore long-read DNA (MapOnt)
/// - PacBio long-read DNA (MapPb)
/// - Illumina CLR long-read DNA (MapIclr)
/// - Long-read RNA with splicing (Splice)
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Preset {
    MapOnt,
    MapPb,
    MapIclr,
    Splice,
}

impl Preset {
    /// Get the preset string used by minimap2
    pub fn as_str(&self) -> &'static str {
        match self {
            Preset::MapOnt => "map-ont",
            Preset::MapPb => "map-pb",
            Preset::MapIclr => "map-iclr",
            Preset::Splice => "splice",
        }
    }
    
    /// Check if this preset is for RNA alignment
    pub fn is_rna(&self) -> bool {
        matches!(self, Preset::Splice)
    }
    
    /// Check if this preset is for DNA alignment
    pub fn is_dna(&self) -> bool {
        matches!(self, Preset::MapOnt | Preset::MapPb | Preset::MapIclr)
    }
}

/// Options for configuring minimap2 alignment
#[derive(Debug, Clone)]
pub struct Minimap2Opts {
    pub genome_dir: PathBuf,
    pub preset: Option<String>,
    pub forward_only: bool, // consider transcript forward strand only ( minimap2 -uf option)
    pub sam_output: bool, // Output in SAM format
    pub extra_params: Vec<String>, // TODO: process additional minimap2 parameters
}

impl Minimap2Opts {
    /// Create new options with the specified genome directory
    pub fn new<P: AsRef<Path>>(genome_dir: P) -> Self {
        Self {
            genome_dir: genome_dir.as_ref().to_path_buf(),
            preset: None,
            forward_only: false,
            sam_output: true, // Default to SAM output
            extra_params: Vec::new(),
        }
    }
    
    /// Set the alignment preset
    pub fn with_preset(mut self, preset: Preset) -> Self {
        self.preset = Some(preset.as_str().to_string());
        self
    }
    
    /// Set forward strand only
    pub fn with_forward_only(mut self, forward_only: bool) -> Self {
        self.forward_only = forward_only;
        self
    }
    
    /// Add extra parameters
    pub fn with_extra_params(mut self, params: Vec<String>) -> Self {
        self.extra_params = params;
        self
    }
    
    // Convert options to C-style arguments for minimap2
    // fn to_c_args(&self) -> Vec<CString> {
    //     let mut args = vec![
    //         CString::new("minimap2").unwrap(),
    //     ];
        
    //     // Add preset if specified
    //     if let Some(ref preset) = self.preset {
    //         args.push(CString::new("-x").unwrap());
    //         args.push(CString::new(preset.as_str()).unwrap());
    //     }
        
    //     // Add SAM output format
    //     if self.sam_output {
    //         args.push(CString::new("-a").unwrap());
    //     }
        
    //     // Add forward only for RNA splice
    //     if self.forward_only && self.preset.as_ref().map_or(false, |p| p == "splice") {
    //         args.push(CString::new("-uf").unwrap());
    //     }
        
    //     // Add extra parameters
    //     for param in &self.extra_params {
    //         args.push(CString::new(param.as_str()).unwrap());
    //     }
        
    //     // Add genome directory
    //     args.push(CString::new(self.genome_dir.to_string_lossy().as_bytes()).unwrap());
        
    //     args
    // }
}

/// A single part of a minimap2 index
#[derive(Debug)]
struct IndexPart {
    inner: *mut minimap2_sys::mm_idx_t,
}

impl Drop for IndexPart {
    fn drop(&mut self) {
        unsafe {
            if !self.inner.is_null() {
                minimap2_sys::mm_idx_destroy(self.inner);
                self.inner = std::ptr::null_mut();
            }
        }
    }
}

/// Minimap2 index that can handle multi-part indices
/// 
/// This structure follows the pattern from example.c with minimap2 where large indices
/// may be split into multiple parts that need to be processed sequentially.
#[derive(Debug)]
pub struct Minimap2Index {
    parts: Vec<IndexPart>, // All index parts loaded from file/sequences
    pub header: Header,
}

impl Minimap2Index {
    /// Load index from the specified options
    /// 
    /// This function can handle both:
    /// - Pre-built index files (.mmi)
    /// - Sequence files (.fasta/.fastq) - will build index on-the-fly
    pub fn load(opts: &Minimap2Opts) -> Result<Self> {
        let path = &opts.genome_dir;
        
        // Check if path exists
        if !path.exists() {
            bail!("Path does not exist: {}", path.display());
        }
        
        let path_cstr = CString::new(path.to_string_lossy().as_bytes())?;
        
        // Initialize index options and mapping options
        let mut iopt = unsafe {
            std::mem::zeroed::<minimap2_sys::mm_idxopt_t>()
        };
        let mut mopt = unsafe {
            std::mem::zeroed::<minimap2_sys::mm_mapopt_t>()
        };
        
        unsafe {
            // First call mm_set_opt with NULL to initialize default values
            minimap2_sys::mm_set_opt(std::ptr::null(), &mut iopt, &mut mopt);
            
            // Then apply preset if specified
            if let Some(ref preset) = opts.preset {
                let preset_cstr = CString::new(preset.as_str())?;
                let result = minimap2_sys::mm_set_opt(preset_cstr.as_ptr(), &mut iopt, &mut mopt);
                if result < 0 {
                    bail!("Failed to apply preset '{}': unknown preset", preset);
                }
            }
        }
        
        // mm_idx_reader_open automatically detects file type and handles both
        // .mmi index files and .fasta/.fastq sequence files
        //
        // type mm_idx_reader_t is the struct for processing provided fasta or index reference file
        let reader = unsafe {
            minimap2_sys::mm_idx_reader_open(
                path_cstr.as_ptr(),
                &mut iopt,
                std::ptr::null(), //path to save index
            )
        };
        
        if reader.is_null() {
            bail!("Failed to open file: {}", path.display());
        }
        
        // Load all parts of the index following example.c pattern:
        // while ((mi = mm_idx_reader_read(r, n_threads)) != 0)
        let mut parts = Vec::new();
        
        loop {
            // type mm_idx_t is designed to contain the true index part
            let idx = unsafe {
                minimap2_sys::mm_idx_reader_read(reader, 4) // to get a part of the index
            };
            
            if idx.is_null() {
                break; // No more parts to read
            }
            
            
            parts.push(IndexPart { inner: idx });
        }
        
        // Close the reader since we've loaded all parts
        unsafe { minimap2_sys::mm_idx_reader_close(reader); }
        
        if parts.is_empty() {
            bail!("Failed to read/build index from: {}", path.display());
        }
        
        // Extract header information from all loaded parts
        let header = Self::extract_header(&parts)?;
        
        Ok(Minimap2Index {
            parts,
            header,
        })
    }
    
    /// Extract header information from all index parts
    /// Already deduplicated and sorted by sequence name
    fn extract_header(parts: &[IndexPart]) -> Result<Header> {
        // Automatically sort by sequence name and handle deduplication using BTreeMap
        let mut sequences = BTreeMap::new();
        
        // Collect all reference sequences from all index parts
        for part in parts {
            let idx = part.inner;
            if idx.is_null() {
                continue;
            }
            
            unsafe {
                let idx_ref = &*idx;
                let n_seq = idx_ref.n_seq;
                
                if idx_ref.seq.is_null() {
                    continue;
                }
                
                // Extract sequence information from this index part
                for i in 0..n_seq {
                    let seq_ptr = idx_ref.seq.add(i as usize);
                    let seq_ref = &*seq_ptr;
                    
                    if seq_ref.name.is_null() || seq_ref.len == 0 {
                        continue;
                    }
                    
                    // Convert C string to Rust string
                    let name_cstr = std::ffi::CStr::from_ptr(seq_ref.name);
                    if let Ok(name) = name_cstr.to_str() {
                        let length = seq_ref.len as usize;
                        // Handle deduplication: if sequence already exists, keep the longer one
                        // or the first one if lengths are equal
                        match sequences.get(name) {
                            Some(&existing_length) => {
                                if length > existing_length {
                                    sequences.insert(name.to_string(), length);
                                }
                            }
                            None => {
                                sequences.insert(name.to_string(), length);
                            }
                        }
                    }
                }
            }
        }
        
        // Build header with deduplicated and sorted sequences
        let mut builder = Header::builder();
        for (name, length) in sequences {
            let length = NonZeroUsize::try_from(length)?;
            let reference_sequence = Map::<ReferenceSequence>::new(length);
            builder = builder.add_reference_sequence(name.as_bytes(), reference_sequence);
        }
        
        Ok(builder.build())
    }
    
    /// Get the number of index parts
    pub fn num_parts(&self) -> usize {
        self.parts.len()
    }
}

unsafe impl Send for Minimap2Index {}
unsafe impl Sync for Minimap2Index {}

/// Inner aligner for implementing thread-safe alignment operations
#[derive(Debug)]
pub struct InnerAligner {
    inner: *mut minimap2_sys::mm_tbuf_t, // thread buffer
    opts: *mut minimap2_sys::mm_mapopt_t, // mapping options
    read_buf: Vec<u8>, // buffer for loading single read sequence
}

impl InnerAligner {
    fn new(opts: &Minimap2Opts) -> Result<Self> {
        let thread_buf = unsafe {
            minimap2_sys::mm_tbuf_init()
        };
        
        if thread_buf.is_null() {
            bail!("Failed to initialize minimap2 thread buffer");
        }
        
        // Set common mapping options, and index-specific parameters will be updated in align_read()
        let map_opts = unsafe {
            // Allocate memory for mapping options
            let opts_ptr = libc::malloc(std::mem::size_of::<minimap2_sys::mm_mapopt_t>()) 
                as *mut minimap2_sys::mm_mapopt_t;
            if opts_ptr.is_null() {
                minimap2_sys::mm_tbuf_destroy(thread_buf);
                bail!("Failed to allocate memory for mapping options");
            }
            
            minimap2_sys::mm_mapopt_init(opts_ptr);
            
            // Apply preset if specified using mm_set_opt
            if let Some(ref preset) = opts.preset {
                Self::apply_preset(opts_ptr, preset)?;
            }
            
            // Set forward only flag for RNA (minimap2 -uf option)
            if opts.forward_only {
                let opts_ref = &mut *opts_ptr;
                // Set MM_F_FOR_ONLY flag (force forward strand only)
                opts_ref.flag |= minimap2_sys::MM_F_FOR_ONLY as i64;
            }
            
            // Enable CIGAR output for SAM format
            if opts.sam_output {
                let opts_ref = &mut *opts_ptr;
                opts_ref.flag |= minimap2_sys::MM_F_CIGAR as i64; // set CIGAR flag in C code
            }
            
            opts_ptr
        };
        
        Ok(InnerAligner {
            inner: thread_buf,
            opts: map_opts,
            read_buf: Vec::with_capacity(1024),
        })
    }
    
    /// Apply preset settings using mm_set_opt
    unsafe fn apply_preset(opts_ptr: *mut minimap2_sys::mm_mapopt_t, preset: &str) -> Result<()> {
        let preset_cstr = CString::new(preset)?;
        
        // Initialize index options for mm_set_opt
        let mut iopt = std::mem::zeroed::<minimap2_sys::mm_idxopt_t>();
        
        // First call mm_set_opt with NULL to initialize default values
        minimap2_sys::mm_set_opt(std::ptr::null(), &mut iopt, opts_ptr);
        
        // Then apply preset using mm_set_opt (example.c pattern)
        let result = minimap2_sys::mm_set_opt(
            preset_cstr.as_ptr(),
            &mut iopt,
            opts_ptr,
        );
        
        if result < 0 {
            bail!("Failed to apply preset '{}': unknown preset", preset);
        }
        
        Ok(())
    }
    
    /// Align a single read. Note that the aligner is mutably borrowed.
    /// We did this on purpose to ensure this function cannot be called concurrently.
    /// Concurrency is implemented in Minimap2Aligner
    /// 
    /// This function follows the pattern from example.c where each part of a 
    /// multi-part index needs to be processed separately.
    fn align_read(&mut self, index: &Minimap2Index, fq: &fastq::Record) -> Result<Vec<RecordBuf>> {
        let query_name = fq.definition().name();
        let sequence = fq.sequence();
        
        // Clear and prepare read buffer
        self.read_buf.clear();
        self.read_buf.extend_from_slice(sequence.as_ref());
        
        let name_cstr = CString::new(query_name.as_ref() as &[u8])?;
        let seq_ptr = self.read_buf.as_ptr() as *const c_char;
        let seq_len = self.read_buf.len() as c_int;
        
        let mut all_alignments = Vec::new();
        
        // Process each index part following example.c pattern
        for part in &index.parts {
            // Update mapping options for this index part (following example.c)
            unsafe {
                minimap2_sys::mm_mapopt_update(self.opts, part.inner); // set index-specific parameters for mapping options
            }

            let mut n_reg = 0i32;
            
            // Call minimap2's core alignment function for this part
            let regs = unsafe {
                minimap2_sys::mm_map(
                    part.inner, // index part
                    seq_len,
                    seq_ptr,
                    &mut n_reg,
                    self.inner, // thread buffer
                    self.opts,
                    name_cstr.as_ptr(),
                )
            };
            
            if !regs.is_null() && n_reg > 0 {
                unsafe {
                    let reg_slice = std::slice::from_raw_parts(regs, n_reg as usize);

                    for reg in reg_slice {
                        if let Ok(record) = self.reg_to_sam_record(&index.header, 
                            std::str::from_utf8(query_name.as_ref())?, 
                            sequence.as_ref(), 
                            reg) {
                            all_alignments.push(record);
                        }
                    }
                    // Free the alignment results for this part
                    libc::free(regs as *mut libc::c_void);
                }
            }
        }
        Ok(all_alignments)
    }
    
    
    /// Convert minimap2 alignment result (mm_reg1_t type) to SAM record
    fn reg_to_sam_record(
        &mut self,
        header: &Header,
        query_name: &str,
        sequence: &[u8],
        reg: &minimap2_sys::mm_reg1_t,
    ) -> Result<RecordBuf> {
        let mut record = RecordBuf::default();
        
        // Set query name
        *record.name_mut() = Some(query_name.as_bytes().to_vec().into());
        
        // Set sequence
        *record.sequence_mut() = sequence.to_vec().into();
        
        // Set reference sequence ID by mapping reg.rid to header reference sequences
        if reg.rid >= 0 {
            let rid = reg.rid as usize;
            // Validate that rid is within the range of reference sequences in header
            if rid < header.reference_sequences().len() {
                *record.reference_sequence_id_mut() = Some(rid);
            }
        }
        
        // Set alignment position (1-based in SAM)
        if reg.rs >= 0 {
            *record.alignment_start_mut() = Some(noodles_core::Position::try_from((reg.rs + 1) as usize)?);
        }
        
        // Set mapping quality
        *record.mapping_quality_mut() = Some(noodles_sam::alignment::record::MappingQuality::try_from(reg.mapq() as u8)?);
        
        // Set SAM flags
        let mut flags = noodles_sam::alignment::record::Flags::empty();
        if reg.rev() != 0 {
            flags |= noodles_sam::alignment::record::Flags::REVERSE_COMPLEMENTED;
        }
        *record.flags_mut() = flags;
        
        // Generate CIGAR operations directly from minimap2 result
        if !reg.p.is_null() {
            if let Ok(cigar_ops) = self.generate_cigar_ops(reg) {
                *record.cigar_mut() = Cigar::from(cigar_ops);
            }
        }
        
        Ok(record)
    }
    
    /// Generate CIGAR operations directly from minimap2 alignment result
    fn generate_cigar_ops(&self, reg: &minimap2_sys::mm_reg1_t) -> Result<Vec<Op>> {
        if reg.p.is_null() {
            return Ok(Vec::new());
        }
        
        let mut ops = Vec::new();
        
        unsafe {
            let p_ref = &*reg.p;
            for i in 0..p_ref.n_cigar {
                // Access CIGAR array safely using as_ptr() and add()
                let cigar_ptr = p_ref.cigar.as_ptr().add(i as usize);
                let cigar_op = *cigar_ptr;
                let op_len = (cigar_op >> 4) as u32;
                let op_type = (cigar_op & 0xf) as u8;
                
                // Map minimap2 CIGAR operation types to noodles_sam Kind
                let kind = match op_type {
                    0 => Kind::Match,           // M - Match/Mismatch
                    1 => Kind::Insertion,       // I - Insertion
                    2 => Kind::Deletion,        // D - Deletion
                    3 => Kind::Skip,            // N - Skip (intron)
                    4 => Kind::SoftClip,        // S - Soft clip
                    5 => Kind::HardClip,        // H - Hard clip
                    7 => Kind::SequenceMatch,   // = - Sequence match
                    8 => Kind::SequenceMismatch, // X - Sequence mismatch
                    _ => Kind::Match,           // Default to match for unknown types
                };
                
                ops.push(Op::new(kind, op_len as usize));
            }
        }
        
        Ok(ops)
    }
}

impl Drop for InnerAligner {
    fn drop(&mut self) {
        unsafe {
            if !self.inner.is_null() {
                minimap2_sys::mm_tbuf_destroy(self.inner);
                self.inner = std::ptr::null_mut();
            }
            if !self.opts.is_null() {
                libc::free(self.opts as *mut libc::c_void);
                self.opts = std::ptr::null_mut();
            }
        }
    }
}

/// Configuration for parallel processing
#[derive(Debug, Clone)]
pub struct ParallelConfig {
    /// Number of threads to use (0 = auto-detect)
    pub num_threads: usize,
    /// Batch size for processing
    pub batch_size: usize,
    /// Enable verbose output
    pub verbose: bool,
}

impl Default for ParallelConfig {
    fn default() -> Self {
        Self {
            num_threads: 0,
            batch_size: 1000,
            verbose: false,
        }
    }
}

/// High-level batch processing interface
/// Minimap2Aligner wrapper with parallel processing configuration
pub struct BatchAligner {
    aligner: Minimap2Aligner,
    config: ParallelConfig,
}

impl BatchAligner {

    pub fn new(opts: Minimap2Opts, config: ParallelConfig) -> Result<Self> {
        let aligner = Minimap2Aligner::new(opts)?;
        Ok(Self { aligner, config })
    }
    
    /// Process a batch of sequences with optimal parallelism
    pub fn process_batch(&self, records: &[fastq::Record]) -> Result<Vec<Vec<RecordBuf>>> {
        self.aligner.align_batch(records, self.config.num_threads)
    }
    
    /// Process sequences from an iterator with streaming
    pub fn process_stream<I, F>(&self, records: I, mut callback: F) -> Result<()>
    where
        I: Iterator<Item = fastq::Record>,
        F: FnMut(Vec<RecordBuf>) -> Result<()>,
    {
        let mut batch = Vec::with_capacity(self.config.batch_size);
        
        for record in records {
            batch.push(record);
            
            if batch.len() >= self.config.batch_size {
                let results = self.process_batch(&batch)?;
                for result in results {
                    callback(result)?;
                }
                batch.clear();
            }
        }
        
        // Process remaining records
        if !batch.is_empty() {
            let results = self.process_batch(&batch)?;
            for result in results {
                callback(result)?;
            }
        }
        Ok(())
    }
    
    /// Get the underlying aligner for direct access
    pub fn aligner(&self) -> &Minimap2Aligner {
        &self.aligner
    }
    
    /// Get configuration
    pub fn config(&self) -> &ParallelConfig {
        &self.config
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    /// Test using .mmi format reference file to initialize Minimap2Aligner
    #[test]
    fn test_load_prebuilt_index() {
        let index_path = "/data/xurui/project/minimap2-rs/data/minimap2_ref/GRCh38.mmi";
        
        let opts = Minimap2Opts::new(index_path)
            .with_preset(Preset::MapOnt);
            
        let result = Minimap2Aligner::new(opts);
        
        match result {
            Ok(aligner) => {
                println!("Successfully loaded prebuilt index from: {}", index_path);

                assert!(aligner.index.parts.len() > 0);
                println!("Number of index parts: {}", aligner.index.parts.len());

                assert!(!aligner.get_header().reference_sequences().is_empty());
                println!("Reference sequences: {:?}", aligner.get_header().reference_sequences());
            }
            Err(e) => {
                println!("Expected failure with fake path: {}", e);
                assert!(e.to_string().contains("does not exist") || 
                        e.to_string().contains("Failed to open file"));
            }
        }
    }

    /// Test using .fasta reference file to build index and initialize Minimap2Aligner
    /// 
    /// NOTE: it will take several minutes to build index from FASTA file
    #[test]
    fn test_dynamic_index_build() {
        let fasta_path = "/data/Public/genome/GRCh38/GRCh38.primary_assembly.genome.fa.gz";
        
        let opts = Minimap2Opts::new(fasta_path)
            .with_preset(Preset::MapPb);
            
        let result = Minimap2Aligner::new(opts);
        
        match result {
            Ok(aligner) => {
                println!("Successfully built dynamic index from: {}", fasta_path);

                assert!(aligner.index.num_parts() > 0);
                println!("Number of index parts: {}", aligner.index.num_parts());

                let header = aligner.get_header();
                assert!(!header.reference_sequences().is_empty());
                
                // Print reference sequence information for debugging
                for (name, ref_seq) in header.reference_sequences() {
                    println!("Reference sequence: {} (length: {})", 
                            String::from_utf8_lossy(name), 
                            ref_seq.length());
                }
            }
            Err(e) => {
                // Expected to fail because the path is fake
                println!("Expected failure with fake path: {}", e);
                assert!(e.to_string().contains("does not exist") || 
                        e.to_string().contains("Failed to open file"));
            }
        }
    }

    use noodles_sam::io::Writer as SamWriter;
    use noodles_sam::alignment::io::Write;
    use std::fs::File;
    use std::path::Path;
    use std::io::BufReader;

    #[test]
    fn test_alignment() -> Result<()> {

        let ref_path = "/data/Public/genome/GRCh38/GRCh38.primary_assembly.genome.fa.gz";
        let reads_path = "/data/xurui/project/minimap2-rs/data/h38_dna/SRR17666426_20.fastq";
        let output_sam_path = "/data/xurui/project/minimap2-rs/data/h38_dna/SRR17666426_20_rs.sam";

        // set expected results from minimap2 (RNAME, POS, FLAG)
        let expected_results = vec![
            ("chr1", 1, 0),
            ("chr1", 1, 0),
            ("chr2", 1, 0),
            ("*", 0, 4),
            ("*", 0, 4),
            ("chr1", 1, 16),
        ];

        // initialize aligner
        let opts = Minimap2Opts::new(&ref_path)
            .with_preset(Preset::MapOnt);
        let mut aligner = Minimap2Aligner::new(opts)?;
        let header = aligner.get_header().clone(); // clone header for later validation

        // read FASTQ file
        let mut reader = fastq::io::Reader::new(BufReader::new(File::open(&reads_path)?));
        let records: Vec<_> = reader.records().collect::<Result<_, _>>()?;

        let mut all_alignments = Vec::new();
        for record in &records {
            let alignments = aligner.align_read(record)?;
            // minimap2 may return multiple alignments for one read, here we only take the first one (usually the best)
            if let Some(best_alignment) = alignments.into_iter().next() {
                all_alignments.push(best_alignment);
            }
        }
        
        assert_eq!(
            all_alignments.len(),
            expected_results.len(),
            "number of alignment results does not match expected"
        );

        for (i, record_buf) in all_alignments.iter().enumerate() {
            let expected = &expected_results[i];
            
            let rname = if let Some(rid) = record_buf.reference_sequence_id() {
                // convert reference ID back to sequence name from header (e.g., "chr1")
                header.reference_sequences().get_index(rid).map(|(name, _)| std::str::from_utf8(name).unwrap()).unwrap_or("*")
            } else {
                "*"
            };

            let pos = record_buf.alignment_start().map_or(0, |p| p.get());
            let flag = record_buf.flags().bits();

            assert_eq!(rname, expected.0, "Record #{} RNAME does not match", i + 1);
            assert_eq!(pos, expected.1, "Record #{} POS does not match", i + 1);
            assert_eq!(flag, expected.2, "Record #{} FLAG does not match", i + 1);
        }
        
        // write SAM file and validate
        let mut sam_writer = SamWriter::new(File::create(&output_sam_path)?);
        sam_writer.write_header(&header)?;
        for record in &all_alignments {
            sam_writer.write_alignment_record(&header, record)?;
        }
        // ensure the content in writer buffer is written to file
        drop(sam_writer); 

        assert!(Path::new(output_sam_path).exists(), "SAM file not created successfully");
        println!("Test successful, SAM file generated at: {}", output_sam_path);
        
        Ok(())
    }

}