fn main() {
    // Allow overriding the library directory via UNO_LIB_DIR environment variable
    let lib_dir = std::env::var("UNO_LIB_DIR").unwrap_or_else(|_| "/usr/local/lib".to_string());
    println!("cargo:rustc-link-search=native={lib_dir}");
    println!("cargo:rustc-link-lib=dylib=uno");
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=UNO_LIB_DIR");
}
