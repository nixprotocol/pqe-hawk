//! Build script: when `cross-check-reference-c` feature is enabled, compiles
//! the reference C (with any debug-trace patches applied) and generates Rust
//! FFI bindings via bindgen. When the feature is disabled, this is a no-op
//! so users of the default crate never need a C toolchain.

fn main() {
    #[cfg(feature = "cross-check-reference-c")]
    build_reference_c();

    #[cfg(not(feature = "cross-check-reference-c"))]
    {
        // No-op when feature is disabled.
    }
}

#[cfg(feature = "cross-check-reference-c")]
fn build_reference_c() {
    use std::path::PathBuf;
    use std::process::Command;

    let c_dir = PathBuf::from("c-reference/hawk-512");
    let patches_dir = PathBuf::from("c-reference/patches");
    let out_dir = PathBuf::from(std::env::var("OUT_DIR").unwrap());

    // Copy the pristine C source into OUT_DIR, apply patches there so the
    // vendored source in the repo stays unmodified.
    let patched_src = out_dir.join("hawk-512-patched");
    if patched_src.exists() {
        std::fs::remove_dir_all(&patched_src).unwrap();
    }
    std::fs::create_dir_all(&patched_src).unwrap();

    for entry in std::fs::read_dir(&c_dir).unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();
        if path.is_file() {
            let dest = patched_src.join(path.file_name().unwrap());
            std::fs::copy(&path, &dest).unwrap();
        }
    }

    // Apply each .patch file in lexicographic order.
    if patches_dir.exists() {
        let mut patches: Vec<_> = std::fs::read_dir(&patches_dir)
            .unwrap()
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.extension().map(|e| e == "patch").unwrap_or(false))
            .collect();
        patches.sort();

        for patch in &patches {
            let status = Command::new("patch")
                .args(["-p1", "-d"])
                .arg(&patched_src)
                .arg("-i")
                .arg(std::fs::canonicalize(patch).unwrap())
                .status()
                .expect("patch command failed (is `patch` installed?)");
            assert!(status.success(), "failed to apply patch: {:?}", patch);
        }
    }

    // Compile all .c files EXCEPT the KAT generator (has its own main()) into
    // a static lib.
    let mut build = cc::Build::new();
    build.include(&patched_src);
    build.flag("-O2");
    build.flag("-DHAWK_DEBUG_TRACE=1");
    // Platform-specific: suppress warnings that shouldn't fail our build.
    build.warnings(false);

    for entry in std::fs::read_dir(&patched_src).unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();
        if path.extension().map(|e| e == "c").unwrap_or(false) {
            let name = path.file_name().unwrap().to_str().unwrap();
            // Exclude standalone programs (have main()).
            if name == "PQCgenKAT_sign.c" {
                continue;
            }
            build.file(&path);
        }
    }
    build.compile("hawk_reference");

    // Generate FFI bindings from the top-level hawk.h header.
    let bindings = bindgen::Builder::default()
        .header(patched_src.join("hawk.h").to_str().unwrap())
        .clang_arg(format!("-I{}", patched_src.display()))
        .clang_arg("-DHAWK_DEBUG_TRACE=1")
        .allowlist_function("hawk_.*")
        .allowlist_function("test_.*")
        .allowlist_type("hawk_.*")
        .generate()
        .expect("bindgen failed");
    bindings
        .write_to_file(out_dir.join("ffi_bindings.rs"))
        .expect("couldn't write bindings");

    println!("cargo:rerun-if-changed=c-reference/");
}
