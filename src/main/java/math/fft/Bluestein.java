/* 
 * Copyright (c) 2017, 2018 Project Nayuki (MIT License) and Stefan Zobel (Apache 2.0)
 * https://www.nayuki.io/page/free-small-fft-in-multiple-languages
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software.
 * - The Software is provided "as is", without warranty of any kind, express or
 *   implied, including but not limited to the warranties of merchantability,
 *   fitness for a particular purpose and noninfringement. In no event shall the
 *   authors or copyright holders be liable for any claim, damages or other
 *   liability, whether in an action of contract, tort or otherwise, arising from,
 *   out of or in connection with the Software or the use or other dealings in the
 *   Software.
 */
/*
 * Any changes, bugfixes or additions made by the maintainers
 * of the https://github.com/stefan-zobel/FFT library are
 * licensed under the Apache License, Version 2.0, as explained
 * at http://www.apache.org/licenses/LICENSE-2.0
 */
package math.fft;

/**
 * Bluestein chirp-z transform
 */
final class Bluestein {

    static ComplexArray forwardDFT(double[] data, double[] imag) {
        // find a power of 2 convolution length m such that m >= n * 2 + 1
        int n = data.length;
        if (n >= 0x20000000) {
            throw new IllegalArgumentException("array too large: " + n);
        }

        double[] cos = new double[n];
        double[] sin = new double[n];
        for (int i = 0; i < n; ++i) {
            int j = (int) ((long) i * i % (n * 2));
            double angle = Math.PI * j / n;
            cos[i] = Math.cos(angle);
            sin[i] = Math.sin(angle);
        }

        int m = Integer.highestOneBit(n) * 4;

        // temporary arrays
        double[] a_re = new double[m];
        double[] a_im = new double[m];
        double[] b_re = new double[m];
        double[] b_im = new double[m];

        b_re[0] = cos[0];
        b_im[0] = sin[0];

        for (int i = 0; i < n; ++i) {
            double sin_i = sin[i];
            double cos_i = cos[i];
            double re_i = data[i];
            double im_i = (imag != null) ? imag[i] : 0.0;
            a_re[i] = re_i * cos_i + im_i * sin_i;
            a_im[i] = -re_i * sin_i + im_i * cos_i;
            if (i != 0) {
                b_re[i] = b_re[m - i] = cos_i;
                b_im[i] = b_im[m - i] = sin_i;
            }
        }

        // convolution
        ComplexArray conv = convolve(new ComplexArray(a_re, a_im, false), new ComplexArray(b_re, b_im, false));
        double[] c_re = conv.re();
        double[] c_im = conv.im();

        // result
        double[] re = new double[n];
        double[] im = new double[n];

        // postprocessing
        for (int i = 0; i < n; ++i) {
            double sin_i = sin[i];
            double cos_i = cos[i];
            double c_re_i = c_re[i];
            double c_im_i = c_im[i];
            double re_i = c_re_i * cos_i + c_im_i * sin_i;
            double im_i = -c_re_i * sin_i + c_im_i * cos_i;
            re[i] = (Math.abs(re_i) <= ComplexArray.TOL) ? 0.0 : re_i;
            im[i] = (Math.abs(im_i) <= ComplexArray.TOL) ? 0.0 : im_i;
        }

        return new ComplexArray(re, im, false);
    }

    static ComplexArray inverseDFT(ComplexArray freqs) {
        ComplexArray inv = forwardDFT(freqs.re(), freqs.im());
        double[] re = inv.re();
        double[] im = inv.im();
        final int n = re.length;
        for (int i = 0; i < n; ++i) {
            double re_i = re[i] / n;
            double im_i = im[i] / n;
            re[i] = (Math.abs(re_i) <= ComplexArray.TOL) ? 0.0 : re_i;
            im[i] = (Math.abs(im_i) <= ComplexArray.TOL) ? 0.0 : im_i;
        }
        for (int i = 1; i <= n / 2; ++i) {
            double re_tmp = re[n - i];
            double im_tmp = im[n - i];
            re[n - i] = re[i];
            re[i] = re_tmp;
            im[n - i] = im[i];
            im[i] = im_tmp;
        }
        return inv;
    }

    private static ComplexArray convolve(ComplexArray x, ComplexArray y) {

        x = Fourier.forwardDFT(x.re(), x.im());
        y = Fourier.forwardDFT(y.re(), y.im());

        double[] x_re = x.re();
        double[] x_im = x.im();
        double[] y_re = y.re();
        double[] y_im = y.im();

        for (int i = 0; i < x_re.length; ++i) {
            double x_re_i = x_re[i];
            double y_re_i = y_re[i];
            double x_im_i = x_im[i];
            double y_im_i = y_im[i];
            x_re[i] = x_re_i * y_re_i - x_im_i * y_im_i;
            x_im[i] = x_im_i * y_re_i + x_re_i * y_im_i;
        }

        return Fourier.inverseDFT(new ComplexArray(x_re, x_im, false));
    }

    private Bluestein() {
        throw new AssertionError();
    }
}
