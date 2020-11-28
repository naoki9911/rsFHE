#include "spqlios-fft.h"
#include <stdio.h>

extern "C" {
    typedef struct {
       FFT_Processor_Spqlios impl;
    } SpqliosImpl;

    SpqliosImpl* Spqlios_new(const int32_t N) {
        FFT_Processor_Spqlios *spqlios = new FFT_Processor_Spqlios(N);
        return (SpqliosImpl *)spqlios;
    }

    void Spqlios_destructor(SpqliosImpl *si) {
        si->impl.~FFT_Processor_Spqlios();
    }

    void Spqlios_ifft_lv1(SpqliosImpl *si, double *res, const uint32_t *src) {
        si->impl.execute_reverse_torus32(res, src);
    }

    void Spqlios_fft_lv1(SpqliosImpl *si, uint32_t *res, const double *src) {
        si->impl.execute_direct_torus32(res, src);
    }
}