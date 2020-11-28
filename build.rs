extern crate cc;
extern crate bindgen;

use std::env;
use std::path::PathBuf;

fn main() {
    println!("cargo:rustc-link-lib=spqlios");

    cc::Build::new()
        .cpp(true)
        .file("src/fft_processor_spqlios.cpp")
        .file("src/spqlios-wrapper.cpp")
        .flag("-std=c++17")
        .flag("-lm")
        .include("src")
        .compile("libspqlios.a");
}