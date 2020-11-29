use std::os::raw::c_double;
use std::os::raw::c_int;
use std::os::raw::c_uint;

pub enum SpqliosImpl {}

extern "C" {
    pub fn Spqlios_new(N: c_int) -> *mut SpqliosImpl;
    pub fn Spqlios_destructor(spqlios: *mut SpqliosImpl);
    pub fn Spqlios_ifft_lv1(spqlios: *mut SpqliosImpl, res: *mut c_double, src: *const c_uint);
    pub fn Spqlios_fft_lv1(spqlios: *mut SpqliosImpl, res: *mut c_uint, src: *const c_double);
}
