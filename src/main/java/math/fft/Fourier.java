/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/*
 * Any changes, bugfixes or additions made by the maintainers
 * of the https://github.com/stefan-zobel/FFT library are
 * licensed under the Apache License, Version 2.0, as explained
 * at http://www.apache.org/licenses/LICENSE-2.0
 */
package math.fft;

/**
 * Derived from org.apache.commons.math3.transform.FastFourierTransformer
 */
public final class Fourier {

    public static ComplexArray forwardDFT(double[] data) {
        final int N = data.length;
        if (N == 0) {
            return new ComplexArray(new double[] {}, new double[] {}, false);
        }
        if (N == 1) {
            return new ComplexArray(data, new double[1], true);
        }
        if (N == 2) {
            double[] dataR = data.clone();
            double srcR0 = dataR[0];
            // X_0 = x_0 + x_1
            dataR[0] = srcR0 + dataR[1];
            // X_1 = x_0 - x_1
            dataR[1] = srcR0 - dataR[1];
            return new ComplexArray(dataR, new double[2], false);
        }
        if (!isPowerOfTwo(N)) {
            return Bluestein.forwardDFT(data, null);
        }
        double[] dataR = data.clone();
        bitReversalShuffle(dataR, null);
        double[] dataI = new double[N];
        fourTermForward(dataR, dataI, N);
        combineEvenOdd(dataR, dataI, N, false);
        postProcess(dataR, dataI, N, false);
        return new ComplexArray(dataR, dataI, false);
    }

    static ComplexArray forwardDFT(double[] real, double[] imag) {
        final int N = real.length;
        if (N == 0) {
            return new ComplexArray(new double[] {}, new double[] {}, false);
        }
        if (N == 1) {
            return new ComplexArray(real, imag, true);
        }
        if (N == 2) {
            double[] dataR = real.clone();
            double[] dataI = imag.clone();
            double srcR0 = dataR[0];
            double srcI0 = dataI[0];
            dataR[0] = srcR0 + dataR[1];
            dataR[1] = srcR0 - dataR[1];
            dataI[0] = srcI0 + dataI[1];
            dataI[1] = srcI0 - dataI[1];
            return new ComplexArray(dataR, dataI, false);
        }
        double[] dataR = real.clone();
        double[] dataI = imag.clone();
        bitReversalShuffle(dataR, dataI);
        fourTermForward(dataR, dataI, N);
        combineEvenOdd(dataR, dataI, N, false);
        postProcess(dataR, dataI, N, false);
        return new ComplexArray(dataR, dataI, false);
    }

    public static ComplexArray inverseDFT(ComplexArray freqs) {
        final int N = freqs.length();
        if (N <= 1) {
            return new ComplexArray(freqs.re(), freqs.im(), true);
        }
        if (N == 2) {
            double[] dataR = freqs.re().clone();
            double[] dataI = freqs.im().clone();
            double srcR0 = dataR[0];
            double srcI0 = dataI[0];
            double srcR1 = dataR[1];
            double srcI1 = dataI[1];
            double scaleFactor = (1.0 / N);
            // X_0 = x_0 + x_1
            dataR[0] = srcR0 + srcR1;
            dataR[0] *= scaleFactor;
            dataI[0] = srcI0 + srcI1;
            dataI[0] *= scaleFactor;
            // X_1 = x_0 - x_1
            dataR[1] = srcR0 - srcR1;
            dataR[1] *= scaleFactor;
            dataI[1] = srcI0 - srcI1;
            dataI[1] *= scaleFactor;
            return new ComplexArray(dataR, dataI, false);
        }
        if (!isPowerOfTwo(N)) {
            return Bluestein.inverseDFT(freqs);
        }
        double[] dataR = freqs.re().clone();
        double[] dataI = freqs.im().clone();
        bitReversalShuffle(dataR, dataI);
        fourTermInverse(dataR, dataI, N);
        combineEvenOdd(dataR, dataI, N, true);
        postProcess(dataR, dataI, N, true);
        return new ComplexArray(dataR, dataI, false);
    }

    private static void fourTermForward(double[] dataR, double[] dataI, int n) {
        for (int i0 = 0; i0 < n; i0 += 4) {
            int i1 = i0 + 1;
            int i2 = i0 + 2;
            int i3 = i0 + 3;

            double srcR0 = dataR[i0];
            double srcI0 = dataI[i0];
            double srcR1 = dataR[i2];
            double srcI1 = dataI[i2];
            double srcR2 = dataR[i1];
            double srcI2 = dataI[i1];
            double srcR3 = dataR[i3];
            double srcI3 = dataI[i3];

            // 4-term DFT
            // X_0 = x_0 + x_1 + x_2 + x_3
            dataR[i0] = srcR0 + srcR1 + srcR2 + srcR3;
            dataI[i0] = srcI0 + srcI1 + srcI2 + srcI3;
            // X_1 = x_0 - x_2 + j * (x_3 - x_1)
            dataR[i1] = srcR0 - srcR2 + (srcI1 - srcI3);
            dataI[i1] = srcI0 - srcI2 + (srcR3 - srcR1);
            // X_2 = x_0 - x_1 + x_2 - x_3
            dataR[i2] = srcR0 - srcR1 + srcR2 - srcR3;
            dataI[i2] = srcI0 - srcI1 + srcI2 - srcI3;
            // X_3 = x_0 - x_2 + j * (x_1 - x_3)
            dataR[i3] = srcR0 - srcR2 + (srcI3 - srcI1);
            dataI[i3] = srcI0 - srcI2 + (srcR1 - srcR3);
        }
    }

    private static void fourTermInverse(double[] dataR, double[] dataI, int n) {
        for (int i0 = 0; i0 < n; i0 += 4) {
            int i1 = i0 + 1;
            int i2 = i0 + 2;
            int i3 = i0 + 3;

            double srcR0 = dataR[i0];
            double srcI0 = dataI[i0];
            double srcR1 = dataR[i2];
            double srcI1 = dataI[i2];
            double srcR2 = dataR[i1];
            double srcI2 = dataI[i1];
            double srcR3 = dataR[i3];
            double srcI3 = dataI[i3];

            // 4-term DFT
            // X_0 = x_0 + x_1 + x_2 + x_3
            dataR[i0] = srcR0 + srcR1 + srcR2 + srcR3;
            dataI[i0] = srcI0 + srcI1 + srcI2 + srcI3;
            // X_1 = x_0 - x_2 + j * (x_3 - x_1)
            dataR[i1] = srcR0 - srcR2 + (srcI3 - srcI1);
            dataI[i1] = srcI0 - srcI2 + (srcR1 - srcR3);
            // X_2 = x_0 - x_1 + x_2 - x_3
            dataR[i2] = srcR0 - srcR1 + srcR2 - srcR3;
            dataI[i2] = srcI0 - srcI1 + srcI2 - srcI3;
            // X_3 = x_0 - x_2 + j * (x_1 - x_3)
            dataR[i3] = srcR0 - srcR2 + (srcI1 - srcI3);
            dataI[i3] = srcI0 - srcI2 + (srcR3 - srcR1);
        }
    }

    private static void combineEvenOdd(double[] dataR, double[] dataI, int n, boolean doInverse) {
        int lastN0 = 4;
        int lastLogN0 = 2;
        while (lastN0 < n) {
            int n0 = lastN0 << 1;
            int logN0 = lastLogN0 + 1;
            double wSubN0R = W_SUB_N_R[logN0];
            double wSubN0I = doInverse ? -W_SUB_N_I[logN0] : W_SUB_N_I[logN0];

            // Combine even/odd transforms of size lastN0 into a transform of
            // size N0 (lastN0 * 2).
            for (int destEvenStartIndex = 0; destEvenStartIndex < n; destEvenStartIndex += n0) {
                int destOddStartIndex = destEvenStartIndex + lastN0;

                double wSubN0ToRR = 1;
                double wSubN0ToRI = 0;

                for (int r = 0; r < lastN0; r++) {
                    double grR = dataR[destEvenStartIndex + r];
                    double grI = dataI[destEvenStartIndex + r];
                    double hrR = dataR[destOddStartIndex + r];
                    double hrI = dataI[destOddStartIndex + r];

                    // dest[destEvenStartIndex + r] = Gr + WsubN0ToR * Hr
                    dataR[destEvenStartIndex + r] = grR + wSubN0ToRR * hrR - wSubN0ToRI * hrI;
                    dataI[destEvenStartIndex + r] = grI + wSubN0ToRR * hrI + wSubN0ToRI * hrR;
                    // dest[destOddStartIndex + r] = Gr - WsubN0ToR * Hr
                    dataR[destOddStartIndex + r] = grR - (wSubN0ToRR * hrR - wSubN0ToRI * hrI);
                    dataI[destOddStartIndex + r] = grI - (wSubN0ToRR * hrI + wSubN0ToRI * hrR);

                    // WsubN0ToR *= WsubN0R
                    double nextWsubN0ToRR = wSubN0ToRR * wSubN0R - wSubN0ToRI * wSubN0I;
                    double nextWsubN0ToRI = wSubN0ToRR * wSubN0I + wSubN0ToRI * wSubN0R;
                    wSubN0ToRR = nextWsubN0ToRR;
                    wSubN0ToRI = nextWsubN0ToRI;
                }
            }

            lastN0 = n0;
            lastLogN0 = logN0;
        }
    }

    private static void postProcess(double[] dataR, double[] dataI, int n, boolean normalize) {
        double scaleFactor = normalize ? (1.0 / n) : 1.0;
        for (int i = 0; i < n; ++i) {
            double re_i = dataR[i] * scaleFactor;
            dataR[i] = (Math.abs(re_i) <= ComplexArray.TOL) ? 0.0 : re_i;
        }
        for (int i = 0; i < n; ++i) {
            double im_i = dataI[i] * scaleFactor;
            dataI[i] = (Math.abs(im_i) <= ComplexArray.TOL) ? 0.0 : im_i;
        }
    }

    /**
     * Performs identical index bit reversal shuffles on two arrays of identical
     * size. Each element in the array is swapped with another element based on
     * the bit-reversal of the index. For example, in an array with length 16,
     * item at binary index 0011 (decimal 3) would be swapped with the item at
     * binary index 1100 (decimal 12).
     *
     * @param a
     *            the first array to be shuffled
     * @param b
     *            the second array to be shuffled
     */
    private static void bitReversalShuffle(double[] a, double[] b) {
        final int n = a.length;
        final int halfOfN = n >> 1;

        int j = 0;
        for (int i = 0; i < n; i++) {
            if (i < j) {
                // swap indices i & j
                double temp = a[i];
                a[i] = a[j];
                a[j] = temp;

                if (b != null) {
                    temp = b[i];
                    b[i] = b[j];
                    b[j] = temp;
                }
            }

            int k = halfOfN;
            while (k <= j && k > 0) {
                j -= k;
                k >>= 1;
            }
            j += k;
        }
    }

    private static boolean isPowerOfTwo(int n) {
        return (n > 0) && ((n & (n - 1)) == 0);
    }

    /**
     * {@code W_SUB_N_R[i]} is the real part of {@code exp(- 2 * i * pi / n)}:
     * {@code W_SUB_N_R[i] = cos(2 * pi/ n)}, where {@code n = 2^i}.
     */
    //@formatter:off
    private static final double[] W_SUB_N_R =
        {  0x1.0p0, -0x1.0p0, 0x1.1a62633145c07p-54, 0x1.6a09e667f3bcdp-1
        , 0x1.d906bcf328d46p-1, 0x1.f6297cff75cbp-1, 0x1.fd88da3d12526p-1, 0x1.ff621e3796d7ep-1
        , 0x1.ffd886084cd0dp-1, 0x1.fff62169b92dbp-1, 0x1.fffd8858e8a92p-1, 0x1.ffff621621d02p-1
        , 0x1.ffffd88586ee6p-1, 0x1.fffff62161a34p-1, 0x1.fffffd8858675p-1, 0x1.ffffff621619cp-1
        , 0x1.ffffffd885867p-1, 0x1.fffffff62161ap-1, 0x1.fffffffd88586p-1, 0x1.ffffffff62162p-1
        , 0x1.ffffffffd8858p-1, 0x1.fffffffff6216p-1, 0x1.fffffffffd886p-1, 0x1.ffffffffff621p-1
        , 0x1.ffffffffffd88p-1, 0x1.fffffffffff62p-1, 0x1.fffffffffffd9p-1, 0x1.ffffffffffff6p-1
        , 0x1.ffffffffffffep-1, 0x1.fffffffffffffp-1, 0x1.0p0, 0x1.0p0
        , 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0
        , 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0
        , 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0
        , 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0
        , 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0
        , 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0
        , 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0
        , 0x1.0p0, 0x1.0p0, 0x1.0p0 };
    //@formatter:on

    /**
     * {@code W_SUB_N_I[i]} is the imaginary part of
     * {@code exp(- 2 * i * pi / n)}: {@code W_SUB_N_I[i] = -sin(2 * pi/ n)},
     * where {@code n = 2^i}.
     */
    //@formatter:off
    private static final double[] W_SUB_N_I =
        {  0x1.1a62633145c07p-52, -0x1.1a62633145c07p-53, -0x1.0p0, -0x1.6a09e667f3bccp-1
        , -0x1.87de2a6aea963p-2, -0x1.8f8b83c69a60ap-3, -0x1.917a6bc29b42cp-4, -0x1.91f65f10dd814p-5
        , -0x1.92155f7a3667ep-6, -0x1.921d1fcdec784p-7, -0x1.921f0fe670071p-8, -0x1.921f8becca4bap-9
        , -0x1.921faaee6472dp-10, -0x1.921fb2aecb36p-11, -0x1.921fb49ee4ea6p-12, -0x1.921fb51aeb57bp-13
        , -0x1.921fb539ecf31p-14, -0x1.921fb541ad59ep-15, -0x1.921fb5439d73ap-16, -0x1.921fb544197ap-17
        , -0x1.921fb544387bap-18, -0x1.921fb544403c1p-19, -0x1.921fb544422c2p-20, -0x1.921fb54442a83p-21
        , -0x1.921fb54442c73p-22, -0x1.921fb54442cefp-23, -0x1.921fb54442d0ep-24, -0x1.921fb54442d15p-25
        , -0x1.921fb54442d17p-26, -0x1.921fb54442d18p-27, -0x1.921fb54442d18p-28, -0x1.921fb54442d18p-29
        , -0x1.921fb54442d18p-30, -0x1.921fb54442d18p-31, -0x1.921fb54442d18p-32, -0x1.921fb54442d18p-33
        , -0x1.921fb54442d18p-34, -0x1.921fb54442d18p-35, -0x1.921fb54442d18p-36, -0x1.921fb54442d18p-37
        , -0x1.921fb54442d18p-38, -0x1.921fb54442d18p-39, -0x1.921fb54442d18p-40, -0x1.921fb54442d18p-41
        , -0x1.921fb54442d18p-42, -0x1.921fb54442d18p-43, -0x1.921fb54442d18p-44, -0x1.921fb54442d18p-45
        , -0x1.921fb54442d18p-46, -0x1.921fb54442d18p-47, -0x1.921fb54442d18p-48, -0x1.921fb54442d18p-49
        , -0x1.921fb54442d18p-50, -0x1.921fb54442d18p-51, -0x1.921fb54442d18p-52, -0x1.921fb54442d18p-53
        , -0x1.921fb54442d18p-54, -0x1.921fb54442d18p-55, -0x1.921fb54442d18p-56, -0x1.921fb54442d18p-57
        , -0x1.921fb54442d18p-58, -0x1.921fb54442d18p-59, -0x1.921fb54442d18p-60 };
    //@formatter:on
}
