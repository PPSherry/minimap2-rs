//! High-level, safe Rust wrapper for minimap2
//!
//! This crate provides a safe, idiomatic interface to minimap2 for long-read sequence alignment.
//! It supports both DNA and RNA alignment with various presets optimized for different sequencing
//! technologies.

use std::ffi::{CStr, CString};
use std::path::{Path, PathBuf};

use noodles_fastq as fastq;
use thiserror::Error;

pub use minimap2_sys;

/// Errors that can occur during alignment operations
#[derive(Error, Debug)]
pub enum AlignmentError {
    #[error("Failed to open index file: {path}")]
    IndexOpen { path: PathBuf },
    
    #[error("Alignment operation failed")]
    AlignmentFailed,
    
    #[error("Invalid sequence data")]
    InvalidSequence,
    
    #[error("Memory allocation failed")]
    OutOfMemory,
    
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    
    #[error("UTF-8 conversion error: {0}")]
    Utf8(#[from] std::str::Utf8Error),
}

pub type Result<T> = std::result::Result<T, AlignmentError>;

/// Alignment presets for different data types and sequencing technologies
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Preset {
    /// Oxford Nanopore long-read DNA (map-ont)
    OntDna,
    /// PacBio long-read DNA (map-pb)
    PacBioDna,
    /// PacBio CLR long-read DNA (map-iclr)  
    PacBioClr,
    /// Long-read RNA with splicing (splice)
    RnaSplice,
    /// Short-read mapping (sr)
    ShortRead,
    /// Assembly-to-assembly, ~5% divergence (asm5)
    Assembly5,
    /// Assembly-to-assembly, ~10% divergence (asm10)
    Assembly10,
    /// Assembly-to-assembly, ~20% divergence (asm20)
    Assembly20,
}

impl Preset {
    /// Get the preset string used by minimap2
    pub fn as_str(self) -> &'static str {
        match self {
            Preset::OntDna => "map-ont",
            Preset::PacBioDna => "map-pb",
            Preset::PacBioClr => "map-iclr",
            Preset::RnaSplice => "splice",
            Preset::ShortRead => "sr",
            Preset::Assembly5 => "asm5",
            Preset::Assembly10 => "asm10",
            Preset::Assembly20 => "asm20",
        }
    }
    
    /// Check if this preset is for RNA alignment
    pub fn is_rna(self) -> bool {
        matches!(self, Preset::RnaSplice)
    }
    
    /// Check if this preset is for DNA alignment
    pub fn is_dna(self) -> bool {
        matches!(self, 
            Preset::OntDna | Preset::PacBioDna | Preset::PacBioClr | Preset::ShortRead
        )
    }
}

/// Alignment result for a single query sequence
#[derive(Debug, Clone)]
pub struct Alignment {
    /// Query sequence name
    pub query_name: String,
    /// Query sequence length
    pub query_len: u32,
    /// Query start position (0-based)
    pub query_start: u32,
    /// Query end position (0-based, exclusive)
    pub query_end: u32,
    /// Target sequence name  
    pub target_name: String,
    /// Target sequence length
    pub target_len: u32,
    /// Target start position (0-based)
    pub target_start: u32,
    /// Target end position (0-based, exclusive)
    pub target_end: u32,
    /// Number of matching residues
    pub matches: u32,
    /// Total alignment block length
    pub block_len: u32,
    /// Mapping quality (0-255)
    pub mapq: u8,
    /// Is alignment on reverse strand?
    pub is_reverse: bool,
    /// CIGAR string if available
    pub cigar: Option<String>,
}

/// Reference sequence information
#[derive(Debug, Clone)]
pub struct ReferenceSequence {
    /// Sequence name
    pub name: String,
    /// Sequence length
    pub length: u32,
}

/// Minimap2 index wrapper
pub struct Index {
    inner: *mut minimap2_sys::mm_idx_t,
    sequences: Vec<ReferenceSequence>,
}

impl Index {
    /// Load index from file
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let path_cstr = CString::new(path.to_string_lossy().as_bytes())
            .map_err(|_| AlignmentError::InvalidSequence)?;
        
        // Open file
        let file_mode = CString::new("rb").unwrap();
        let file_ptr = unsafe {
            libc::fopen(path_cstr.as_ptr(), file_mode.as_ptr())
        };
        
        if file_ptr.is_null() {
            return Err(AlignmentError::IndexOpen { 
                path: path.to_path_buf() 
            });
        }
        
        // Load index
        let idx = unsafe {
            let result = minimap2_sys::mm_idx_load(file_ptr as *mut minimap2_sys::FILE);
            libc::fclose(file_ptr);
            result
        };
        
        if idx.is_null() {
            return Err(AlignmentError::IndexOpen { 
                path: path.to_path_buf() 
            });
        }
        
        // Extract reference sequences
        let sequences = Self::extract_sequences(idx)?;
        
        Ok(Index {
            inner: idx,
            sequences,
        })
    }
    
    /// Get reference sequences
    pub fn sequences(&self) -> &[ReferenceSequence] {
        &self.sequences
    }
    
    fn extract_sequences(idx: *mut minimap2_sys::mm_idx_t) -> Result<Vec<ReferenceSequence>> {
        let mut sequences = Vec::new();
        
        unsafe {
            let idx_ref = &*idx;
            for i in 0..idx_ref.n_seq {
                let seq_ptr = *idx_ref.seq.offset(i as isize);
                let name = CStr::from_ptr(seq_ptr.name).to_str()?.to_string();
                let length = seq_ptr.len;
                sequences.push(ReferenceSequence { name, length });
            }
        }
        
        Ok(sequences)
    }
}

impl Drop for Index {
    fn drop(&mut self) {
        if !self.inner.is_null() {
            unsafe {
                minimap2_sys::mm_idx_destroy(self.inner);
            }
        }
    }
}

unsafe impl Send for Index {}
unsafe impl Sync for Index {}

/// Minimap2 aligner
pub struct Aligner {
    index: Index,
    opts: *mut minimap2_sys::mm_mapopt_t,
    thread_buffer: *mut minimap2_sys::mm_tbuf_t,
}

impl Aligner {
    /// Create a new aligner with the specified preset and index
    pub fn new<P: AsRef<Path>>(index_path: P, preset: Preset) -> Result<Self> {
        Self::builder().preset(preset).build(index_path)
    }
    
    /// Create an aligner builder for advanced configuration
    pub fn builder() -> AlignerBuilder {
        AlignerBuilder::new()
    }
    
    /// Get reference sequences from the index
    pub fn sequences(&self) -> &[ReferenceSequence] {
        self.index.sequences()
    }
    
    /// Align a single sequence
    pub fn align(&mut self, record: &fastq::Record) -> Result<Vec<Alignment>> {
        let query_name = record.definition().name().to_string();
        let sequence = record.sequence();
        
        // Convert to C strings
        let name_cstr = CString::new(query_name.as_bytes())
            .map_err(|_| AlignmentError::InvalidSequence)?;
        let seq_cstr = CString::new(sequence.as_ref() as &[u8])
            .map_err(|_| AlignmentError::InvalidSequence)?;
        
        let mut n_reg = 0i32;
        let regs = unsafe {
            minimap2_sys::mm_map(
                self.index.inner,
                sequence.len() as i32,
                seq_cstr.as_ptr(),
                &mut n_reg,
                self.thread_buffer,
                self.opts,
                name_cstr.as_ptr(),
            )
        };
        
        if regs.is_null() && n_reg > 0 {
            return Err(AlignmentError::AlignmentFailed);
        }
        
        // Convert results
        let mut alignments = Vec::new();
        if n_reg > 0 {
            let reg_slice = unsafe {
                std::slice::from_raw_parts(regs, n_reg as usize)
            };
            
            for reg in reg_slice {
                if let Ok(alignment) = self.reg_to_alignment(&query_name, sequence.len(), reg) {
                    alignments.push(alignment);
                }
            }
        }
        
        // Free memory
        if !regs.is_null() {
            unsafe {
                libc::free(regs as *mut libc::c_void);
            }
        }
        
        Ok(alignments)
    }
    
    /// Align paired-end sequences
    pub fn align_pair(
        &mut self, 
        r1: &fastq::Record, 
        r2: &fastq::Record
    ) -> Result<(Vec<Alignment>, Vec<Alignment>)> {
        let alignments1 = self.align(r1)?;
        let alignments2 = self.align(r2)?;
        Ok((alignments1, alignments2))
    }
    
    fn reg_to_alignment(
        &self,
        query_name: &str,
        query_len: usize,
        reg: &minimap2_sys::mm_reg1_t,
    ) -> Result<Alignment> {
        let target_name = if reg.rid >= 0 && (reg.rid as usize) < self.index.sequences.len() {
            self.index.sequences[reg.rid as usize].name.clone()
        } else {
            "*".to_string()
        };
        
        let target_len = if reg.rid >= 0 && (reg.rid as usize) < self.index.sequences.len() {
            self.index.sequences[reg.rid as usize].length
        } else {
            0
        };
        
        let cigar = if !reg.p.is_null() {
            self.extract_cigar(reg)?
        } else {
            None
        };
        
        Ok(Alignment {
            query_name: query_name.to_string(),
            query_len: query_len as u32,
            query_start: reg.qs as u32,
            query_end: reg.qe as u32,
            target_name,
            target_len,
            target_start: reg.rs as u32,
            target_end: reg.re as u32,
            matches: reg.mlen as u32,
            block_len: reg.blen as u32,
            mapq: reg.mapq() as u8,
            is_reverse: reg.rev() != 0,
            cigar,
        })
    }
    
    fn extract_cigar(&self, reg: &minimap2_sys::mm_reg1_t) -> Result<Option<String>> {
        if reg.p.is_null() {
            return Ok(None);
        }
        
        let mut cigar = String::new();
        unsafe {
            let p_ref = &*reg.p;
            for i in 0..p_ref.n_cigar {
                let cigar_ptr = p_ref.cigar.as_ptr().add(i as usize);
                let cigar_op = *cigar_ptr;
                let op_len = cigar_op >> 4;
                let op_type = (cigar_op & 0xf) as u8;
                
                let op_char = match op_type {
                    0 => 'M', // Match
                    1 => 'I', // Insertion
                    2 => 'D', // Deletion
                    3 => 'N', // Skip
                    4 => 'S', // Soft clip
                    5 => 'H', // Hard clip
                    7 => '=', // Sequence match
                    8 => 'X', // Sequence mismatch
                    _ => 'M', // Default
                };
                
                cigar.push_str(&format!("{}{}", op_len, op_char));
            }
        }
        
        Ok(Some(cigar))
    }
}

impl Drop for Aligner {
    fn drop(&mut self) {
        unsafe {
            if !self.thread_buffer.is_null() {
                minimap2_sys::mm_tbuf_destroy(self.thread_buffer);
            }
            if !self.opts.is_null() {
                libc::free(self.opts as *mut libc::c_void);
            }
        }
    }
}

unsafe impl Send for Aligner {}

/// Builder for configuring aligner options
pub struct AlignerBuilder {
    preset: Preset,
    min_chain_score: Option<i32>,
    min_dp_score: Option<i32>,
    bandwidth: Option<i32>,
    max_gap: Option<i32>,
    threads: Option<i32>,
    forward_only: bool,
    output_cigar: bool,
}

impl AlignerBuilder {
    fn new() -> Self {
        Self {
            preset: Preset::OntDna,
            min_chain_score: None,
            min_dp_score: None,
            bandwidth: None,
            max_gap: None,
            threads: None,
            forward_only: false,
            output_cigar: true,
        }
    }
    
    /// Set alignment preset
    pub fn preset(mut self, preset: Preset) -> Self {
        self.preset = preset;
        self
    }
    
    /// Set minimum chaining score
    pub fn min_chain_score(mut self, score: i32) -> Self {
        self.min_chain_score = Some(score);
        self
    }
    
    /// Set minimum DP alignment score
    pub fn min_dp_score(mut self, score: i32) -> Self {
        self.min_dp_score = Some(score);
        self
    }
    
    /// Set alignment bandwidth
    pub fn bandwidth(mut self, bw: i32) -> Self {
        self.bandwidth = Some(bw);
        self
    }
    
    /// Set maximum gap length
    pub fn max_gap(mut self, gap: i32) -> Self {
        self.max_gap = Some(gap);
        self
    }
    
    /// Set number of threads (0 for auto)
    pub fn threads(mut self, n: i32) -> Self {
        self.threads = Some(n);
        self
    }
    
    /// Force forward strand only (useful for RNA)
    pub fn forward_only(mut self, forward: bool) -> Self {
        self.forward_only = forward;
        self
    }
    
    /// Enable/disable CIGAR output
    pub fn output_cigar(mut self, cigar: bool) -> Self {
        self.output_cigar = cigar;
        self
    }
    
    /// Build the aligner with the specified index
    pub fn build<P: AsRef<Path>>(self, index_path: P) -> Result<Aligner> {
        // Load index
        let index = Index::from_file(index_path)?;
        
        // Allocate and configure mapping options
        let mm_opts = unsafe {
            let size = std::mem::size_of::<minimap2_sys::mm_mapopt_t>();
            let opts_ptr = libc::malloc(size) as *mut minimap2_sys::mm_mapopt_t;
            if opts_ptr.is_null() {
                return Err(AlignmentError::OutOfMemory);
            }
            
            // Initialize with defaults
            minimap2_sys::mm_mapopt_init(opts_ptr);
            
            // Apply preset
            minimap2_sys::mm_mapopt_update(opts_ptr, index.inner);
            
            // Apply custom settings
            let opts_ref = &mut *opts_ptr;
            
            if let Some(min_chain_score) = self.min_chain_score {
                opts_ref.min_chain_score = min_chain_score;
            }
            if let Some(min_dp_score) = self.min_dp_score {
                opts_ref.min_dp_max = min_dp_score;
            }
            if let Some(bandwidth) = self.bandwidth {
                opts_ref.bw = bandwidth;
            }
            if let Some(max_gap) = self.max_gap {
                opts_ref.max_gap = max_gap;
            }
            
            // Set flags
            if self.preset.is_rna() {
                opts_ref.flag |= minimap2_sys::MM_F_SPLICE as i64;
            }
            if self.forward_only {
                opts_ref.flag |= minimap2_sys::MM_F_FOR_ONLY as i64;
            }
            if self.output_cigar {
                opts_ref.flag |= minimap2_sys::MM_F_CIGAR as i64;
            }
            
            opts_ptr
        };
        
        // Initialize thread buffer
        let thread_buffer = unsafe {
            minimap2_sys::mm_tbuf_init()
        };
        
        if thread_buffer.is_null() {
            unsafe {
                libc::free(mm_opts as *mut libc::c_void);
            }
            return Err(AlignmentError::OutOfMemory);
        }
        
        Ok(Aligner {
            index,
            opts: mm_opts,
            thread_buffer,
        })
    }
}

impl Default for AlignerBuilder {
    fn default() -> Self {
        Self::new()
    }
}