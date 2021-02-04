/* 
 * Bluestein chirp-z transform test
 * 
 * Copyright (c) 2017 Project Nayuki. (MIT License)
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

import static org.junit.Assert.*;

import org.junit.Test;

import math.fft.ComplexArray;
import math.fft.Fourier;

public class BluesteinForwardOddTest {

    private static final int NUM_TESTS = 50;

    @Test
    public void testForward() {
        double maxLogErr = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < NUM_TESTS; ++i) {
            int size = TestUtils.randLengthOdd();
            System.out.println("Testing odd forward " + size + " ...");
            double[] data = TestUtils.randomData(size);
            double logErr = testForward(data);
            maxLogErr = Math.max(logErr, maxLogErr);
        }
        System.out.printf("\nMax log err = %.1f%n", maxLogErr);
        System.out.println("Test " + (maxLogErr < -10.0 ? "passed" : "failed"));
        assertTrue(maxLogErr < -10.0);
    }

    private static double testForward(double[] data) {
        long start = System.currentTimeMillis();
        ComplexArray result = Fourier.forwardDFT(data);
        long end = System.currentTimeMillis();

        long start2 = System.currentTimeMillis();
        ComplexArray expected = ComplexArray.naiveForwarDFT(data);
        long end2 = System.currentTimeMillis();

        double err = TestUtils.log10RmsError(expected, result);
        System.out.print(data.length + " took: " + (end - start) + " ms vs. " + (end2 - start2) + " ms | logerr = ");
        System.out.printf("%5.1f%n\n", err);
        return err;
    }
}
