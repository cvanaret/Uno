fn main() {
    println!("cargo:rustc-link-search=native=/usr/local/lib");
    println!("cargo:rustc-link-lib=dylib=uno");
    // ensure cargo re-runs this if the library changes
    println!("cargo:rerun-if-changed=build.rs");
}
