/*
 * Copyright 2018 Stefan Zobel
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package math.fft;

/**
 * For computations with arrays of complex numbers.
 * <p>
 * Note that indexes in {@link #set(int, double, double)} are 1-based!
 */
public final class ComplexArray {

    /** The IEEE 754 machine epsilon from Cephes: {@code (2^-53)} */
    private static final double MACH_EPS = 1.11022302462515654042e-16;
    static final double TOL = 5.0 * MACH_EPS;
    private static final double TWO_PI = 2.0 * Math.PI;

    private final double[] re;
    private final double[] im;

    public ComplexArray(int size) {
        if (size < 0) {
            throw new IllegalArgumentException("size < 0 : " + size);
        }
        re = new double[size];
        im = new double[size];
    }

    public ComplexArray(double[] re) {
        this(re.clone(), new double[re.length], false);
    }

    public ComplexArray(double[] re, double[] im) {
        this(re, im, true);
    }

    public ComplexArray(double[] re, double[] im, boolean copy) {
        if (re.length != im.length) {
            throw new IllegalArgumentException(re.length + " != " + im.length);
        }
        if (copy) {
            this.re = re.clone();
            this.im = im.clone();
        } else {
            this.re = re;
            this.im = im;
        }
    }

    public void set(int index, double re, double im) {
        checkArg(index);
        this.re[index - 1] = re;
        this.im[index - 1] = im;
    }

    public ComplexArray naiveForwardDFT() {
        return naiveDFT(-1.0, re, im, 1.0);
    }

    public ComplexArray naiveInverseDFT() {
        return naiveDFT(1.0, re, im, (1.0 / re.length));
    }

    private static ComplexArray naiveDFT(double sign, double[] re, double[] im, double scale) {
        int N = re.length;
        double[] imag = new double[N];
        double[] real = new double[N];
        double[] cos = new double[N];
        double[] sin = new double[N];
        for (int i = 0; i < N; ++i) {
            double angle = (sign * TWO_PI * i) / N;
            cos[i] = Math.cos(angle);
            sin[i] = Math.sin(angle);
        }
        for (int i = 0; i < N; ++i) {
            double rZ = 0.0;
            double iZ = 0.0;
            for (long j = 0; j < N; ++j) {
                int idx = (int) ((i * j) % N);
                double cosine = cos[idx];
                double sine = sin[idx];
                double rX = re[(int) j];
                double iY = (im == null) ? 0.0 : im[(int) j];
                rZ += cosine * rX - sine * iY;
                iZ += sine * rX + cosine * iY;
            }
            double x = scale * rZ;
            double y = scale * iZ;
            if (Math.abs(x) <= TOL) {
                x = 0.0;
            }
            if (Math.abs(y) <= TOL) {
                y = 0.0;
            }
            real[i] = x;
            imag[i] = y;
        }
        return new ComplexArray(real, imag, false);
    }

    public static ComplexArray naiveForwarDFT(double[] data) {
        return naiveDFT(-1.0, data, null, 1.0);
    }

    public static ComplexArray naiveInverseDFT(ComplexArray freqs) {
        return naiveDFT(1.0, freqs.re, freqs.im, (1.0 / freqs.re.length));
    }

    public double[] absSquared() {
        return absSquaredScaled(false);
    }

    // for power density spectrum
    public double[] absSquaredScaled() {
        return absSquaredScaled(true);
    }

    private double[] absSquaredScaled(boolean withScaling) {
        double[] real = re;
        double[] imag = im;
        int N = real.length;
        double[] res = new double[N];
        double scale = withScaling ? N : 1.0;
        for (int i = 0; i < N; ++i) {
            double rX = real[i];
            double iY = imag[i];
            double square = (rX * rX + iY * iY) / scale;
            if (square <= TOL) {
                square = 0.0;
            }
            res[i] = square;
        }
        return res;
    }

    public ComplexArray fftshift() {
        return shift(false);
    }

    public ComplexArray ifftshift() {
        return shift(true);
    }

    private ComplexArray shift(boolean inverse) {
        final int length = re.length;
        int mid = -1;
        double[] re_this = re;
        double[] im_this = im;
        double[] re_shift = new double[length];
        double[] im_shift = new double[length];
        if (length % 2 == 0) {
            mid = (length / 2);
            System.arraycopy(re_this, 0, re_shift, mid, mid);
            System.arraycopy(re_this, mid, re_shift, 0, mid);
            System.arraycopy(im_this, 0, im_shift, mid, mid);
            System.arraycopy(im_this, mid, im_shift, 0, mid);
        } else {
            mid = (length - 1) / 2;
            if (inverse) {
                System.arraycopy(re_this, 0, re_shift, mid + 1, mid);
                System.arraycopy(re_this, mid, re_shift, 0, mid + 1);
                System.arraycopy(im_this, 0, im_shift, mid + 1, mid);
                System.arraycopy(im_this, mid, im_shift, 0, mid + 1);
            } else {
                System.arraycopy(re_this, 0, re_shift, mid, mid + 1);
                System.arraycopy(re_this, mid + 1, re_shift, 0, mid);
                System.arraycopy(im_this, 0, im_shift, mid, mid + 1);
                System.arraycopy(im_this, mid + 1, im_shift, 0, mid);
            }
        }
        return new ComplexArray(re_shift, im_shift, false);
    }

    public static double[] dot(ComplexArray a, ComplexArray b) {
        if (a.length() != b.length()) {
            throw new IllegalArgumentException("Unequal dimensions: " + a.length() + " != " + b.length());
        }
        if (a.length() == 0) {
            throw new IllegalArgumentException("Arrays are empty: length = 0");
        }
        double res_re = 0.0;
        double res_im = 0.0;
        double[] a_re_ = a.re;
        double[] b_re_ = b.re;
        double[] a_im_ = a.im;
        double[] b_im_ = b.im;
        for (int i = 0; i < a_re_.length; ++i) {
            double a_re = a_re_[i];
            double b_re = b_re_[i];
            double a_im = a_im_[i];
            double b_im = b_im_[i];
            double re_i = a_re * b_re - a_im * b_im;
            double im_i = a_re * b_im + a_im * b_re;
            re_i = (Math.abs(re_i) <= TOL) ? 0.0 : re_i;
            im_i = (Math.abs(im_i) <= TOL) ? 0.0 : im_i;
            res_re += re_i;
            res_im += im_i;
        }
        res_re = (Math.abs(res_re) <= TOL) ? 0.0 : res_re;
        res_im = (Math.abs(res_im) <= TOL) ? 0.0 : res_im;
        return new double[] { res_re, res_im };
    }

    public static ComplexArray elementwiseProduct(ComplexArray a, ComplexArray b) {
        if (a.length() != b.length()) {
            throw new IllegalArgumentException("Unequal dimensions: " + a.length() + " != " + b.length());
        }
        double[] real = new double[a.length()];
        double[] imag = new double[a.length()];
        double[] a_re_ = a.re;
        double[] b_re_ = b.re;
        double[] a_im_ = a.im;
        double[] b_im_ = b.im;
        for (int i = 0; i < a_re_.length; ++i) {
            double a_re = a_re_[i];
            double b_re = b_re_[i];
            double a_im = a_im_[i];
            double b_im = b_im_[i];
            double re_i = a_re * b_re - a_im * b_im;
            double im_i = a_re * b_im + a_im * b_re;
            re_i = (Math.abs(re_i) <= TOL) ? 0.0 : re_i;
            im_i = (Math.abs(im_i) <= TOL) ? 0.0 : im_i;
            real[i] = re_i;
            imag[i] = im_i;
        }
        return new ComplexArray(real, imag, false);
    }

    public double[] re() {
        return re;
    }

    public double[] im() {
        return im;
    }

    public int length() {
        return re.length;
    }

    public String toString() {
        int max = length() - 1;
        if (max == -1) {
            return "[]";
        }
        StringBuilder b = new StringBuilder(40 * (max + 1));
        b.append('[');
        for (int i = 0; ; i++) {
            b.append(re[i]).append("  ").append(im[i]).append('i');
            if (i == max) {
                return b.append(']').toString();
            }
            b.append(",\n ");
        }
    }

    private void checkArg(int idx) {
        if (idx < 1 || idx > re.length) {
            throw new IllegalArgumentException("Invalid index " + idx + " for [1.." + re.length + "] array");
        }
    }
}
