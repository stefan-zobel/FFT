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

import static org.junit.Assert.*;

import org.junit.Test;

import math.fft.ComplexArray;

/**
 * Test for {@link ComplexArray#fftshift()} and {@link ComplexArray#ifftshift()}
 */
public class ComplexArrayFFTShiftTest {

    @Test
    public void testOddSize() throws Exception {
        int SIZE = 7;
        ComplexArray ca1 = new ComplexArray(SIZE);
        for (int i = 0; i < SIZE; ++i) {
            ca1.set(i + 1, i + 1.0, -(i + 1.0));
        }
        System.out.println(ca1.toString());
        System.out.println("----------\n");
        double[] re1 = ca1.re();
        double[] im1 = ca1.im();
        assertArrayEquals(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 }, re1, 0.0);
        assertArrayEquals(new double[] { -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0 }, im1, 0.0);

        ComplexArray ca2 = ca1.fftshift();
        System.out.println(ca2.toString());
        System.out.println("----------\n");
        double[] re2 = ca2.re();
        double[] im2 = ca2.im();
        assertArrayEquals(new double[] { 5.0, 6.0, 7.0, 1.0, 2.0, 3.0, 4.0 }, re2, 0.0);
        assertArrayEquals(new double[] { -5.0, -6.0, -7.0, -1.0, -2.0, -3.0, -4.0 }, im2, 0.0);

        ComplexArray ca3 = ca2.ifftshift();
        System.out.println(ca3.toString());
        System.out.println("----------\n");
        double[] re3 = ca3.re();
        double[] im3 = ca3.im();
        assertArrayEquals(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 }, re3, 0.0);
        assertArrayEquals(new double[] { -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0 }, im3, 0.0);
        System.out.println("testOddSize success\n");
    }

    @Test
    public void testEvenSize() throws Exception {
        int SIZE = 6;
        ComplexArray ca1 = new ComplexArray(SIZE);
        for (int i = 0; i < SIZE; ++i) {
            ca1.set(i + 1, i + 1.0, -(i + 1.0));
        }
        System.out.println(ca1.toString());
        System.out.println("----------\n");
        double[] re1 = ca1.re();
        double[] im1 = ca1.im();
        assertArrayEquals(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 }, re1, 0.0);
        assertArrayEquals(new double[] { -1.0, -2.0, -3.0, -4.0, -5.0, -6.0 }, im1, 0.0);

        ComplexArray ca2 = ca1.fftshift();
        System.out.println(ca2.toString());
        System.out.println("----------\n");
        double[] re2 = ca2.re();
        double[] im2 = ca2.im();
        assertArrayEquals(new double[] { 4.0, 5.0, 6.0, 1.0, 2.0, 3.0 }, re2, 0.0);
        assertArrayEquals(new double[] { -4.0, -5.0, -6.0, -1.0, -2.0, -3.0 }, im2, 0.0);

        ComplexArray ca3 = ca2.ifftshift();
        System.out.println(ca3.toString());
        System.out.println("----------\n");
        double[] re3 = ca3.re();
        double[] im3 = ca3.im();
        assertArrayEquals(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 }, re3, 0.0);
        assertArrayEquals(new double[] { -1.0, -2.0, -3.0, -4.0, -5.0, -6.0 }, im3, 0.0);
        System.out.println("testEvenSize success\n");
    }
}
