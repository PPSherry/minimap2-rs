
use std::env;
use std::path::PathBuf;

fn main() {
    // 1. Set C source code path
    let minimap2_src_dir = PathBuf::from("minimap2");
    println!("cargo:rerun-if-changed={}", minimap2_src_dir.display());

    // 2. Define base object files based on Makefile OBJS
    let base_sources = vec![
        "kthread.c", "kalloc.c", "misc.c", "bseq.c", "sketch.c", "sdust.c",
        "options.c", "index.c", "lchain.c", "align.c", "hit.c", "seed.c",
        "jump.c", "map.c", "format.c", "pe.c", "esterr.c", "splitidx.c",
    ];

    let c_sources: Vec<PathBuf> = base_sources
        .iter()
        .map(|&name| minimap2_src_dir.join(name))
        .collect();

    // 3. Detect target architecture for SIMD support
    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap();
    let is_aarch64 = target_arch == "aarch64";
    let is_arm = target_arch == "arm";
    let is_x86_64 = target_arch == "x86_64";
    let is_x86 = target_arch == "x86";
    
    let arm_neon = is_aarch64 || is_arm;
    let sse2only = env::var("CARGO_FEATURE_SSE2ONLY").is_ok();

    // 4. Compile base sources
    let mut builder = cc::Build::new();
    builder
        .files(&c_sources)
        .include(&minimap2_src_dir)
        .flag("-Wall")
        .flag("-O2")
        .flag("-Wc++-compat")
        .define("HAVE_KALLOC", None)
        .define("cputime", "mm_cputime")
        .define("realtime", "mm_realtime")
        .warnings(false);

    // 5. Add architecture-specific compiler flags and compile SIMD files separately
    if arm_neon {
        builder.include(minimap2_src_dir.join("sse2neon"));
        if !is_aarch64 {
            builder.flag("-D_FILE_OFFSET_BITS=64");
            builder.flag("-mfpu=neon");
            builder.flag("-fsigned-char");
        } else {
            builder.flag("-D_FILE_OFFSET_BITS=64");
            builder.flag("-fsigned-char");
        }
        
        // Compile NEON versions (using SSE source files with NEON defines)
        let mut neon_builder = cc::Build::new();
        neon_builder
            .include(&minimap2_src_dir)
            .include(minimap2_src_dir.join("sse2neon"))
            .flag("-Wall")
            .flag("-O2")
            .flag("-Wc++-compat")
            .define("HAVE_KALLOC", None)
            .define("KSW_SSE2_ONLY", None)
            .define("__SSE2__", None)
            .warnings(false);
            
        if !is_aarch64 {
            neon_builder.flag("-D_FILE_OFFSET_BITS=64");
            neon_builder.flag("-mfpu=neon");
            neon_builder.flag("-fsigned-char");
        } else {
            neon_builder.flag("-D_FILE_OFFSET_BITS=64");
            neon_builder.flag("-fsigned-char");
        }
        
        neon_builder
            .files(vec![
                minimap2_src_dir.join("ksw2_extz2_sse.c"),
                minimap2_src_dir.join("ksw2_extd2_sse.c"),
                minimap2_src_dir.join("ksw2_exts2_sse.c"),
            ])
            .compile("minimap2_neon");
        
    } else {
        // x86/x86_64 SSE support
        if is_x86_64 || is_x86 {
            builder.flag_if_supported("-msse2");
        }
        
        // Compile ksw2_ll_sse.c with SSE2
        let mut sse_ll_builder = cc::Build::new();
        sse_ll_builder
            .include(&minimap2_src_dir)
            .flag("-Wall")
            .flag("-O2")
            .flag("-Wc++-compat")
            .define("HAVE_KALLOC", None)
            .flag_if_supported("-msse2")
            .warnings(false)
            .file(minimap2_src_dir.join("ksw2_ll_sse.c"))
            .compile("minimap2_sse_ll");
        
        if sse2only {
            // Compile SSE2-only versions
            let mut sse2_builder = cc::Build::new();
            sse2_builder
                .include(&minimap2_src_dir)
                .flag("-Wall")
                .flag("-O2")
                .flag("-Wc++-compat")
                .define("HAVE_KALLOC", None)
                .flag_if_supported("-msse2")
                .warnings(false)
                .files(vec![
                    minimap2_src_dir.join("ksw2_extz2_sse.c"),
                    minimap2_src_dir.join("ksw2_extd2_sse.c"),
                    minimap2_src_dir.join("ksw2_exts2_sse.c"),
                ])
                .compile("minimap2_sse2");
        } else {
            // Compile multiple SSE versions for CPU dispatch
            
            // SSE4.1 versions
            let mut sse41_builder = cc::Build::new();
            sse41_builder
                .include(&minimap2_src_dir)
                .flag("-Wall")
                .flag("-O2")
                .flag("-Wc++-compat")
                .define("HAVE_KALLOC", None)
                .define("KSW_CPU_DISPATCH", None)
                .flag_if_supported("-msse4.1")
                .warnings(false)
                .files(vec![
                    minimap2_src_dir.join("ksw2_extz2_sse.c"),
                    minimap2_src_dir.join("ksw2_extd2_sse.c"),
                    minimap2_src_dir.join("ksw2_exts2_sse.c"),
                ])
                .compile("minimap2_sse41");
            
            // SSE2 fallback versions (same source files, different flags and output names)
            let mut sse2_fallback_builder = cc::Build::new();
            sse2_fallback_builder
                .include(&minimap2_src_dir)
                .flag("-Wall")
                .flag("-O2")
                .flag("-Wc++-compat")
                .define("HAVE_KALLOC", None)
                .define("KSW_CPU_DISPATCH", None)
                .define("KSW_SSE2_ONLY", None)
                .flag_if_supported("-msse2")
                .flag_if_supported("-mno-sse4.1")
                .warnings(false)
                .files(vec![
                    minimap2_src_dir.join("ksw2_extz2_sse.c"),
                    minimap2_src_dir.join("ksw2_extd2_sse.c"),
                    minimap2_src_dir.join("ksw2_exts2_sse.c"),
                ])
                .compile("minimap2_sse2_fallback");
            
            // Dispatch mechanism
            let mut dispatch_builder = cc::Build::new();
            dispatch_builder
                .include(&minimap2_src_dir)
                .flag("-Wall")
                .flag("-O2")
                .flag("-Wc++-compat")
                .define("HAVE_KALLOC", None)
                .define("KSW_CPU_DISPATCH", None)
                .flag_if_supported("-msse4.1")
                .warnings(false)
                .file(minimap2_src_dir.join("ksw2_dispatch.c"))
                .compile("minimap2_dispatch");
        }
    }

    // 8. Link with required libraries  
    println!("cargo:rustc-link-lib=z");
    println!("cargo:rustc-link-lib=m");
    if cfg!(target_os = "linux") || cfg!(target_os = "macos") {
        println!("cargo:rustc-link-lib=pthread");
    }

    // 9. Compile the main library
    builder.compile("minimap2");

    // 10. Configure and run bindgen to generate Rust bindings
    println!("cargo:rerun-if-changed=minimap2/minimap.h");
    let mut bindgen_builder = bindgen::Builder::default()
        .header("minimap2/minimap.h")
        // Only generate bindings for public APIs to avoid namespace pollution
        .allowlist_function("mm_.*")
        .allowlist_type("mm_.*")
        .allowlist_var("MM_.*")
        // include header directory & define HAVE_KALLOC
        .clang_arg("-Iminimap2")
        .clang_arg("-DHAVE_KALLOC")
        // Preserve documentation comments from C code
        .generate_comments(true)
        // Add useful derives to generated types
        .derive_default(true)
        .derive_debug(true)
        .derive_copy(true);

    // Add architecture-specific clang args
    if arm_neon {
        bindgen_builder = bindgen_builder.clang_arg("-Iminimap2/sse2neon");
        if !is_aarch64 {
            bindgen_builder = bindgen_builder
                .clang_arg("-D_FILE_OFFSET_BITS=64")
                .clang_arg("-mfpu=neon")
                .clang_arg("-fsigned-char");
        } else {
            bindgen_builder = bindgen_builder
                .clang_arg("-D_FILE_OFFSET_BITS=64")
                .clang_arg("-fsigned-char");
        }
    } else if is_x86_64 || is_x86 {
        bindgen_builder = bindgen_builder.clang_arg("-msse2");
        if !sse2only {
            bindgen_builder = bindgen_builder.clang_arg("-msse4.1");
        }
    }

    let bindings = bindgen_builder
        .generate()
        .expect("Unable to generate bindings");

    // 11. Write generated bindings to $OUT_DIR/bindings.rs
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");

}
