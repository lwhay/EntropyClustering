/*
 * Copyright 2007-2016 by The Regents of the Wuhan University of China.
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * you may obtain a copy of the License from
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package ics.whu.edu.cn.madrix.basis;

import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

import ics.whu.edu.cn.madrix.common.exceptions.MadrixException;

public class MadrixUtils {
    public static final double EPS = 1e-7;
    public static final double PI = Math.PI;

    public static int[] vectorDiff(double[] src, int sb, int se, double[] des, int db, int de) throws MadrixException {
        if (se - sb != de - db) {
            throw new MadrixException("Vector diff with different intput");
        }
        int[] diff = new int[se - sb];
        for (int i = 0; i < se - sb; i++) {
            diff[i] = (src[sb + i] == des[db + i]) ? 0 : 1;
        }
        return diff;
    }

    public static double vectorSum(long[] vector) {
        double sum = .0;
        for (int i = 0; i < vector.length; i++) {
            sum += vector[i];
        }
        return sum;
    }

    public static void vectorSet(double[] vector, double set) {
        for (int i = 0; i < vector.length; i++) {
            vector[i] = set;
        }
    }

    public static void matrixSet(double[][] matrix, double set) {
        for (int i = 0; i < matrix.length; i++) {
            vectorSet(matrix[i], set);
        }
    }

    public static double vectorSum(int[] vector) {
        double sum = .0;
        for (int i = 0; i < vector.length; i++) {
            sum += vector[i];
        }
        return sum;
    }

    public static double[] matrixRowSum(double[][] matrix) {
        double[] rsum = new double[matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                rsum[j] += matrix[i][j];
            }
        }
        return rsum;
    }

    public static double[] matrixColSum(double[][] matrix) {
        double[] csum = new double[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            csum[i] = vectorSum(matrix[i]);
        }
        return csum;
    }

    public static double[] matrixColSum(int[][] matrix) {
        double[] csum = new double[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            csum[i] = vectorSum(matrix[i]);
        }
        return csum;
    }

    public static void matrixColNorm(double[][] matrix) {
        double[] csum = matrixColSum(matrix);
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                matrix[i][j] /= csum[i];
            }
        }
    }

    public static void matrixRowNorm(double[][] matrix) {
        double[] rsum = matrixRowSum(matrix);
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                matrix[i][j] /= rsum[j];
            }
        }
    }

    //    public static double[] matrixRowSum(double[][] matrix) {
    //        double[] rsum = new double[matrix[0].length];
    //        for (int i = 0; i < matrix[0].length; i++) {
    //            for (int j = 0; j < matrix.length; j++) {
    //                rsum[i] += matrix[j][i];
    //            }
    //        }
    //        return rsum;
    //    }

    public static int vectorMax(double[] vector) {
        int pos = -1;
        double min = Double.MIN_VALUE;
        for (int i = 0; i < vector.length; i++) {
            if (vector[i] > min) {
                pos = i;
                min = vector[pos];
            }
        }
        return pos;
    }

    public static int vectorMax(int[] vector) {
        int pos = -1;
        int min = Integer.MIN_VALUE;
        for (int i = 0; i < vector.length; i++) {
            if (vector[i] > min) {
                pos = i;
                min = vector[pos];
            }
        }
        return pos;
    }

    public static int vectorMin(double[] vector) {
        int pos = -1;
        double min = Double.MAX_VALUE;
        for (int i = 0; i < vector.length; i++) {
            if (vector[i] < min) {
                pos = i;
                min = vector[pos];
            }
        }
        return pos;
    }

    public static double vectorSum(double[] vector) {
        double sum = .0;
        for (int i = 0; i < vector.length; i++) {
            sum += vector[i];
        }
        return sum;
    }

    public static double vectorMpl(double[] vector) {
        double sum = .0;
        for (int i = 0; i < vector.length; i++) {
            sum *= vector[i];
        }
        return sum;
    }

    public static boolean doubleEqual(double src, double des) {
        return (Math.abs(src - des) < EPS);
    }

    public static boolean doubleEqual(double src, double des, double precision) {
        double eps = Math.pow(0.1, (precision + 1)) * 5;
        return (Math.abs(src - des) < eps);
    }

    public static double[] getColumnVector(List<List<Double>> matrix, int index) throws MadrixException {
        if (matrix.size() < 1 || matrix.get(0).size() < index) {
            throw new MadrixException("Invalid matrix column");
        }
        double[] vector = new double[matrix.size()];
        for (int i = 0; i < matrix.size(); i++) {
            vector[i] = matrix.get(i).get(index);
        }
        return vector;
    }

    public static double[] getColumnVector(double[][] matrix, int index) throws MadrixException {
        if (matrix.length < 1 || matrix[0].length < index) {
            throw new MadrixException("Invalid matrix column");
        }
        double[] vector = new double[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            vector[i] = matrix[i][index];
        }
        return vector;
    }

    public static double[] getSubVector(double[] vector, int start, int length) throws MadrixException {
        if (start + length > vector.length) {
            throw new MadrixException("Error on subVector: " + start + length + " out of: " + vector.length);
        }
        double[] ret = new double[length];
        System.arraycopy(vector, start, ret, 0, length);
        return ret;
    }

    public static double[][] getSubMadrix(double[][] matrix, int sRow, int nRows, int sCol, int nCols)
            throws MadrixException {
        if (sRow + nRows > matrix.length || sCol + nCols > matrix[0].length) {
            throw new MadrixException("Extracting larger than the matirx. " + " rows: " + sRow + nRows + " cols: "
                    + sCol + nCols + " out of: " + matrix.length + " * " + matrix[0].length);
        }
        double[][] ret = new double[nRows][];
        for (int i = 0; i < nRows; i++) {
            ret[i] = new double[nCols];
            System.arraycopy(matrix[i + sRow], sCol, ret, 0, nCols);
        }
        return ret;
    }

    public static double[][] genMatrix(int dim, int len) {
        double[][] matrix = new double[dim][];
        for (int i = 0; i < dim; i++) {
            matrix[i] = new double[len];
        }
        return matrix;
    }

    public static double[][] genMatrixAndInit(int dim, int len, double init) {
        double[][] matrix = new double[dim][];
        for (int i = 0; i < dim; i++) {
            matrix[i] = new double[len];
            Arrays.fill(matrix[i], init);
        }
        return matrix;
    }

    public static double[][] genMatrixWithInnerProduct(double[] vector) {
        double[][] matrix = new double[vector.length][];
        for (int i = 0; i < vector.length; i++) {
            matrix[i] = dupVector(vector);
            mplVector(matrix[i], vector[i]);
        }
        return matrix;
    }

    public static double[] genVectorAndInit(int len, double init) {
        double[] vector = genVector(len);
        Arrays.fill(vector, init);
        return vector;
    }

    public static double[] genUnitVector(int len) {
        double[] vector = genVector(len);
        Arrays.fill(vector, 1);
        return vector;
    }

    public static double[][] genUnitMatrix(int dim, int card) {
        int min = Math.min(dim, card);
        double[][] matrix = genMatrix(dim, card);
        for (int i = 0; i < min; i++) {
            matrix[i][i] = 1;
        }
        return matrix;
    }

    public static double[] genAscVector(int len) {
        double[] vector = new double[len];
        for (int i = 0; i < len; i++) {
            vector[i] = i;
        }
        return vector;
    }

    public static double[] genVector(int len) {
        double[] vector = new double[len];
        return vector;
    }

    public static double[][] dupMatrix(double[][] matrix) {
        double[][] dup = new double[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            dup[i] = dupVector(matrix[i]);
        }
        return dup;
    }

    //For copy the sub-columns from a given matrix.
    //Both row and column will be copied with [begin, end].
    public static void copyMatrix(double[][] dec, double[][] src, int rows, int cols, int rowb, int rowe, int colb,
            int cole) throws MadrixException {
        if (rows + rowe - rowb + 1 > dec.length)
            throw new MadrixException("copyVector overflows: " + "start->" + rows + " length->" + (rowe - rowb + 1)
                    + " actual->" + dec.length);
        int card = rowe - rowb + 1;
        for (int i = 0; i < card; i++) {
            try {
                copyVector(dec[i + rows], src[i + rowb], cols, colb, cole);
            } catch (MadrixException e) {
                throw new MadrixException(e);
            }
        }
    }

    public static void copyVector(double[] dec, double[] src, int start, int begin, int end) throws MadrixException {
        if (start + end - begin + 1 > dec.length)
            throw new MadrixException("copyVector overflows: " + "start->" + start + " length->" + (end - begin + 1)
                    + " actual->" + dec.length);
        System.arraycopy(src, begin, dec, start, end - begin + 1);
    }

    //For duplicating the sub-columns from a given matrix.
    public static double[][] dupMatrix(double[][] matrix, int begin, int end) {
        double[][] ret = new double[matrix.length][];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = dupVector(matrix[i], begin, end);
        }
        return ret;
    }

    //Both row and column will be duplicated with [begin, end].
    public static double[][] dupMatrix(double[][] matrix, int rowb, int rowe, int colb, int cole) {
        int card = rowe - rowb + 1;
        double[][] ret = new double[card][];
        for (int i = 0; i < card; i++) {
            ret[i] = dupVector(matrix[i + rowb], colb, cole);
        }
        return ret;
    }

    public static double[][] dupMatrix(double[][] matrix, int[] perm, int begin, int end) {
        int card = perm.length;
        double[][] ret = new double[card][];
        for (int i = 0; i < card; i++) {
            ret[i] = MadrixUtils.dupVector(matrix[perm[i]], begin, end);
        }
        return ret;
    }

    public static double[] dupVector(double[] vector) {
        return Arrays.copyOf(vector, vector.length);
    }

    public static double[] dupVector(double[] vector, int begin, int end) {
        return Arrays.copyOfRange(vector, begin, end + 1);
    }

    public static double[][] extendVectorToMatrix(double[] vector, int card) {
        double[][] matrix = new double[card][];
        for (int i = 0; i < card; i++) {
            matrix[i] = dupVector(vector);
        }
        return matrix;
    }

    public static double[][] matrixConcate(double[][] des, double[][] src) throws MadrixException {
        if (des[0].length != src[0].length) {
            throw new MadrixException("Dimensions of the two concated matrixes are different.");
        }
        double[][] mat = new double[des.length + src.length][];
        for (int i = 0; i < des.length; i++) {
            mat[i] = des[i];
        }
        for (int i = 0; i < src.length; i++) {
            mat[des.length + i] = src[i];
        }
        return mat;
    }

    public static double[][] matrixExtendVectors(double[][] des, double[] row, double[] col) throws MadrixException {
        if (des.length != col.length || des[0].length != row.length) {
            throw new MadrixException("Extending the madrix with mismatching dims: " + " row * col: " + des.length
                    + " * " + des[0].length + " by row: " + col.length + " col: " + row.length);
        }
        double[][] ret = new double[des.length + 1][des[0].length + 1];
        for (int i = 0; i < des.length; i++) {
            System.arraycopy(des[i], 0, ret[i], 0, des[0].length);
            ret[i][des[i].length] = col[i];
        }
        System.arraycopy(row, 0, ret[des.length], 0, row.length);
        return ret;
    }

    public static double[][] matrixExtend(double[][] des, double[] vector) throws MadrixException {
        if (des.length != vector.length) {
            throw new MadrixException("Dimentsion of the madrix and column are different");
        }
        double[][] mat = new double[des.length][des[0].length + 1];
        for (int i = 0; i < des.length; i++) {
            mat[i] = new double[des[0].length + 1];
            System.arraycopy(des[i], 0, mat[i], 0, des[0].length);
            mat[i][des[i].length] = vector[i];
        }
        return mat;
    }

    public static double[][] matrixAppend(double[][] des, double[] vector) throws MadrixException {
        if (des[0].length != vector.length) {
            throw new MadrixException("Dimensions of the madrix and row are different");
        }
        double[][] mat = new double[des.length + 1][];
        for (int i = 0; i < des.length; i++) {
            mat[i] = des[i];
        }
        mat[mat.length - 1] = vector;
        return mat;
    }

    public static void matrixInsertVector(double[][] des, double[] vector, int idx) throws MadrixException {
        if (des[0].length != vector.length) {
            throw new MadrixException("Dimensions of the madrix and vector are different");
        }
        System.arraycopy(vector, 0, des[idx], 0, vector.length);
    }

    public static double[] vectorConcate(double[] des, double[] src) throws MadrixException {
        if (des.length != src.length) {
            throw new MadrixException("Dimensions of the two concated vectors are different.");
        }
        double[] vector = new double[des.length + src.length];
        System.arraycopy(des, 0, vector, 0, des.length);
        System.arraycopy(src, 0, vector, des.length, src.length);
        return vector;
    }

    public static double[] vectorAppend(double[] des, double src) {
        double[] vector = new double[des.length + 1];
        System.arraycopy(des, 0, vector, 0, des.length);
        vector[vector.length - 1] = src;
        return vector;
    }

    //Simple find
    public static int gtVector(double[] vector, int begin, double value) throws MadrixException {
        if (begin < 0 || begin > vector.length) {
            throw new MadrixException("Find out of range.");
        }
        int at = -1;
        for (int i = 0; i < vector.length - begin; i++) {
            if (vector[begin + i] > value) {
                at = begin + i;
                break;
            }
        }
        return at;
    }

    public static int ltVector(double[] vector, int begin, double value) throws MadrixException {
        if (begin < 0 || begin > vector.length) {
            throw new MadrixException("Find out of range.");
        }
        int at = -1;
        for (int i = 0; i < vector.length - begin; i++) {
            if (vector[begin + i] < value) {
                at = begin + i;
                break;
            }
        }
        return at;
    }

    public static int geVector(double[] vector, int begin, double value) throws MadrixException {
        if (begin < 0 || begin > vector.length) {
            throw new MadrixException("Find out of range.");
        }
        int at = -1;
        for (int i = 0; i < vector.length - begin; i++) {
            if (vector[begin + i] >= value) {
                at = begin + i;
                break;
            }
        }
        return at;
    }

    public static int leVector(double[] vector, int begin, double value) throws MadrixException {
        if (begin < 0 || begin > vector.length) {
            throw new MadrixException("Find out of range.");
        }
        int at = -1;
        for (int i = 0; i < vector.length - begin; i++) {
            if (vector[begin + i] <= value) {
                at = begin + i;
                break;
            }
        }
        return at;
    }

    public static BitSet gtVectorAll(double[] vector, int begin, double value) throws MadrixException {
        BitSet bs = new BitSet();
        if (begin < 0 || begin > vector.length) {
            throw new MadrixException("Find out of range.");
        }
        for (int i = 0; i < vector.length - begin; i++) {
            if (vector[begin + i] > value) {
                bs.set(begin + i);
            }
        }
        return bs;
    }

    public static BitSet[] gtMatrixAll(double[][] vector, double value) throws MadrixException {
        BitSet[] bs = new BitSet[vector.length];
        for (int i = 0; i < vector.length; i++) {
            bs[i] = gtVectorAll(vector[i], 0, value);
        }
        return bs;
    }

    public static BitSet flattenSets(BitSet[] bss) {
        BitSet bs = new BitSet();
        for (int i = 0; i < bss.length; i++) {
            for (int idx = 0; idx >= 0; idx = bs.nextSetBit(idx + 1)) {
                bs.set(i * bss[0].length() + idx);
            }
        }
        return bs;
    }

    public static BitSet ltVectorAll(double[] vector, int begin, double value) throws MadrixException {
        BitSet bs = new BitSet();
        if (begin < 0 || begin > vector.length) {
            throw new MadrixException("Find out of range.");
        }
        for (int i = 0; i < vector.length - begin; i++) {
            if (vector[begin + i] < value) {
                bs.set(begin + i);
            }
        }
        return bs;
    }

    public static BitSet geVectorAll(double[] vector, int begin, double value) throws MadrixException {
        BitSet bs = new BitSet();
        if (begin < 0 || begin > vector.length) {
            throw new MadrixException("Find out of range.");
        }
        for (int i = 0; i < vector.length - begin; i++) {
            if (vector[begin + i] >= value) {
                bs.set(begin + i);
            }
        }
        return bs;
    }

    public static BitSet leVectorAll(double[] vector, int begin, double value) throws MadrixException {
        BitSet bs = new BitSet();
        if (begin < 0 || begin > vector.length) {
            throw new MadrixException("Find out of range.");
        }
        for (int i = 0; i < vector.length - begin; i++) {
            if (vector[begin + i] <= value) {
                bs.set(begin + i);
            }
        }
        return bs;
    }

    //Matrix triangle functions with generation.
    public static double[] sinVectorGenerate(double[] vector) {
        double[] sin = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            sin[i] = Math.sin(vector[i]);
        }
        return sin;
    }

    public static double[] cosVectorGenerate(double[] vector) {
        double[] cos = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            cos[i] = Math.cos(vector[i]);
        }
        return cos;
    }

    public static double[] tanVectorGenerate(double[] vector) {
        double[] tan = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            tan[i] = Math.tan(vector[i]);
        }
        return tan;
    }

    public static double[] cotVectorGenerate(double[] vector) {
        double[] cot = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            cot[i] = 1 / Math.tan(vector[i]);
        }
        return cot;
    }

    public static double[][] sinMatrixGenerate(double[][] matrix) {
        double[][] sin = new double[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            sin[i] = sinVectorGenerate(matrix[i]);
        }
        return sin;
    }

    public static double[][] cosMatrixGenerate(double[][] matrix) {
        double[][] cos = new double[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            cos[i] = cosVectorGenerate(matrix[i]);
        }
        return cos;
    }

    public static double[][] tanMatrixGenerate(double[][] matrix) {
        double[][] tan = new double[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            tan[i] = tanVectorGenerate(matrix[i]);
        }
        return tan;
    }

    public static double[][] cotMatrixGenerate(double[][] matrix) {
        double[][] cot = new double[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            cot[i] = cotVectorGenerate(matrix[i]);
        }
        return cot;
    }

    //Matrix triangle functions with replacement.
    public static void sinVector(double[] vector) {
        for (int i = 0; i < vector.length; i++) {
            vector[i] = Math.sin(vector[i]);
        }
    }

    public static void cosVector(double[] vector) {
        for (int i = 0; i < vector.length; i++) {
            vector[i] = Math.cos(vector[i]);
        }
    }

    public static void tanVector(double[] vector) {
        for (int i = 0; i < vector.length; i++) {
            vector[i] = Math.tan(vector[i]);
        }
    }

    public static void cotVector(double[] vector) {
        for (int i = 0; i < vector.length; i++) {
            vector[i] = 1 / Math.tan(vector[i]);
        }
    }

    public static void sinMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            sinVector(matrix[i]);
        }
    }

    public static void cosMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            cosVector(matrix[i]);
        }
    }

    public static void tanMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            tanVector(matrix[i]);
        }
    }

    public static void cotMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            cotVector(matrix[i]);
        }
    }

    //Matrix polynomial with generation.
    public static double[] powVectorGenerate(double[] vector, double exp) {
        double[] pow = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            pow[i] = Math.pow(vector[i], exp);
        }
        return pow;
    }

    public static double[][] powMatrixGenerate(double[][] matrix, double exp) {
        double[][] pow = new double[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            pow[i] = powVectorGenerate(matrix[i], exp);
        }
        return pow;
    }

    public static double[] expVectorGenerate(double[] vector) {
        double[] exp = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            exp[i] = Math.exp(vector[i]);
        }
        return exp;
    }

    public static double[][] expMatrixGenerate(double[][] matrix) {
        double[][] exp = new double[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            exp[i] = expVectorGenerate(matrix[i]);
        }
        return exp;
    }

    //Matrix polynomial with replacement.
    public static void powVector(double[] vector, double exp) {
        for (int i = 0; i < vector.length; i++) {
            vector[i] = Math.pow(vector[i], exp);
        }
    }

    public static void powMatrix(double[][] matrix, double exp) {
        for (int i = 0; i < matrix.length; i++) {
            powVector(matrix[i], exp);
        }
    }

    public static void expVector(double[] vector) {
        for (int i = 0; i < vector.length; i++) {
            vector[i] = Math.exp(vector[i]);
        }
    }

    public static void powMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            expVector(matrix[i]);
        }
    }

    //Constant matrix algebrick with replacement.
    public static void addVector(double[] vector, double value) {
        for (int i = 0; i < vector.length; i++) {
            vector[i] = vector[i] + value;
        }
    }

    public static void subVector(double[] vector, double value) {
        addVector(vector, -value);
    }

    public static void mplVector(double[] vector, double value) {
        for (int i = 0; i < vector.length; i++) {
            vector[i] *= value;
        }
    }

    public static void divVector(double[] vector, double value) {
        for (int i = 0; i < vector.length; i++) {
            vector[i] /= value;
        }
    }

    public static void absVector(double[] vector) {
        for (int i = 0; i < vector.length; i++) {
            vector[i] = Math.abs(vector[i]);
        }
    }

    public static void addMatrix(double[][] matrix, double value) {
        for (int i = 0; i < matrix.length; i++) {
            addVector(matrix[i], value);
        }
    }

    public static void subMatrix(double[][] matrix, double value) {
        for (int i = 0; i < matrix.length; i++) {
            subVector(matrix[i], value);
        }
    }

    public static void mplMatrix(double[][] matrix, double value) {
        for (int i = 0; i < matrix.length; i++) {
            mplVector(matrix[i], value);
        }
    }

    public static void divMatrix(double[][] matrix, double value) {
        for (int i = 0; i < matrix.length; i++) {
            divVector(matrix[i], value);
        }
    }

    public static void absMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            absVector(matrix[i]);
        }
    }

    //Constant diag polynomial with replacement.
    public static void addMatrixDiag(double[][] matrix, double value) throws MadrixException {
        if (matrix.length != matrix[0].length) {
            throw new MadrixException("Diag on a non-quare matrix.");
        }
        for (int i = 0; i < matrix.length; i++) {
            matrix[i][i] += value;
        }
    }

    public static void subMatrixDiag(double[][] matrix, double value) throws MadrixException {
        addMatrixDiag(matrix, -value);
    }

    public static void mplMatrixDiag(double[][] matrix, double value) throws MadrixException {
        if (matrix.length != matrix[0].length) {
            throw new MadrixException("Diag on a non-quare matrix.");
        }
        for (int i = 0; i < matrix.length; i++) {
            matrix[i][i] *= value;
        }
    }

    public static void divMatrixDiag(double[][] matrix, double value) throws MadrixException {
        if (matrix.length != matrix[0].length) {
            throw new MadrixException("Diag on a non-quare matrix.");
        }
        for (int i = 0; i < matrix.length; i++) {
            matrix[i][i] /= value;
        }
    }

    //Constant matrix algbrick with generation.
    public static double[] subVectorGeneration(double[] vector, double value) {
        double[] sub = Arrays.copyOf(vector, vector.length);
        subVector(sub, value);
        return sub;
    }

    public static double[][] subMatrixGenerator(double[][] matrix, double value) {
        double[][] sub = new double[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            sub[i] = subVectorGeneration(matrix[i], value);
        }
        return sub;
    }

    public static double[] mplVectorGeneration(double[] vector, double value) {
        double[] mpl = Arrays.copyOf(vector, vector.length);
        mplVector(mpl, value);
        return mpl;
    }

    public static double[][] mplMatrixGeneration(double[][] matrix, double value) {
        double[][] mpl = new double[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            mpl[i] = mplVectorGeneration(matrix[i], value);
        }
        return mpl;
    }

    //Matrix operation with replacement
    public static void matrixSub(double[][] dest, double[][] subs) throws MadrixException {
        if (dest.length != subs.length || dest[0].length != subs[0].length) {
            throw new MadrixException("Destination dimension non-equal to the substitutor");
        }
        for (int i = 0; i < dest.length; i++) {
            for (int j = 0; j < dest[i].length; j++) {
                dest[i][j] -= subs[i][j];
            }
        }
    }

    public static void matrixSum(double[][] des, double[][] src) throws MadrixException {
        if (des.length != src.length || des[0].length != src[0].length) {
            throw new MadrixException("Destination dimension non-equal to the summation");
        }
        for (int i = 0; i < des.length; i++) {
            for (int j = 0; j < des[i].length; j++) {
                des[i][j] += src[i][j];
            }
        }
    }

    public static double[][] matrixSumToNew(double[][] des, double[][] src) throws MadrixException {
        if (des.length != src.length || des[0].length != src[0].length) {
            throw new MadrixException("Destination dimension non-equal to the summation");
        }
        double[][] re = new double[des.length][des[0].length];
        for (int i = 0; i < des.length; i++) {
            for (int j = 0; j < des[i].length; j++) {
                re[i][j] = des[i][j] + src[i][j];
            }
        }
        return re;
    }

    public static void matrixSum(double[][] des, double value) {
        for (int i = 0; i < des.length; i++) {
            try {
                vectorSum(des[i], value);
            } catch (MadrixException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
    }

    //General matrix multiply, left division, right division, and beyond.
    public static double[][] matrixMpl(double[][] dest, double[][] src) throws MadrixException {
        if (dest.length != src[0].length || dest[0].length != src.length) {
            throw new MadrixException("Destination dimension does not match the src on multiply.");
        }
        double[][] matrix = genMatrix(dest.length, dest.length);
        for (int i = 0; i < dest.length; i++) {
            for (int j = 0; j < dest.length; j++) {
                for (int k = 0; k < src.length; k++) {
                    matrix[i][j] += dest[i][k] * src[k][j];
                }
            }
        }
        return matrix;
    }

    public static double[][] matrixMpl_col(double[][] dest, double[][] src) throws MadrixException {
        int m = dest.length;
        int n = dest[0].length;
        int p = src[0].length;
        if (dest[0].length != src.length) {
            throw new MadrixException("Destination dimension does not match the src on multiply.");
        }
        double[][] matrix = genMatrix(m, p);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < p; j++) {
                matrix[i][j] = Muladd(dest[i], src, j, n);
            }
        }
        return matrix;
    }

    public static double Muladd(double[] a, double[][] b, int j, int n) {
        double sum = 0;
        for (int k = 0; k < n; k++) {
            sum += a[k] * b[k][j];
        }
        return sum;
    }

    //    public static double[][] TransposeMatrix(double[][] matrix) {
    //        int col = matrix.length;
    //        int row = matrix[0].length;
    //        double[][] re = new double[row][col];
    //        for (int i = 0; i < row; i++) {
    //            for (int j = 0; j < col; j++) {
    //                re[i][j] = matrix[j][i];
    //            }
    //        }
    //        return re;
    //    }

    public static int[][] TransposeMatrix(int[][] matrix) {
        int col = matrix.length;
        int row = matrix[0].length;
        int[][] re = new int[row][col];
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                re[i][j] = matrix[j][i];
            }
        }
        return re;
    }

    public static double[][] rectangleInverse(double[][] src) throws MadrixException {
        double[][] ret = new double[src.length][];
        for (int i = 0; i < src.length; i++) {
            ret[i] = MadrixUtils.genVector(src[i].length);
            for (int j = 0; j < src[i].length; j++) {
                ret[i][j] = (i == j ? 1D : 0D);
            }
        }
        return ret;
    }

    public static double[][] rectangleDiagonalMpl(double[][] dest, double[][] src) throws MadrixException {
        if (dest[0].length != src.length) {
            throw new MadrixException("Destination dimension does not match the src on multiply.");
        }
        double[][] matrix = genMatrix(src.length, src[0].length);
        for (int i = 0; i < src.length; i++) {
            for (int j = 0; j < src[0].length; j++) {
                for (int k = 0; k < src.length; k++) {
                    if (i == k) {
                        matrix[i][j] += dest[i][k] * src[k][j];
                    }
                }
            }
        }
        return matrix;
    }

    public static double[][] rectangleMpl(double[][] dest, double[][] src) throws MadrixException {
        if (dest[0].length != src.length) {
            throw new MadrixException("Destination dimension does not match the src on multiply.");
        }
        double[][] matrix = genMatrix(src.length, src[0].length);
        for (int i = 0; i < src.length; i++) {
            for (int j = 0; j < src[0].length; j++) {
                for (int k = 0; k < src.length; k++) {
                    matrix[i][j] += dest[i][k] * src[k][j];
                }
            }
        }
        return matrix;
    }

    public static void matrixLeftDiv(double[][] dest, double[][] src) throws MadrixException {
        if (dest.length != src.length || dest[0].length != src[0].length) {
            throw new MadrixException("Destination dimensions not-equal for left division");
        }
        matrixDiv(src, dest);
    }

    public static double[][] matrixLeftDivGenerate(double[][] dest, double[][] src) throws MadrixException {
        if (dest.length != src.length || dest[0].length != src[0].length) {
            throw new MadrixException("Destination dimensions not-equal for left division generation");
        }
        double[][] matrix = dupMatrix(src);
        matrixDiv(matrix, dest);
        return matrix;
    }

    public static void matrixRightDiv(double[][] dest, double[][] src) throws MadrixException {
        if (dest.length != src.length || dest[0].length != src[0].length) {
            throw new MadrixException("Destination dimensions not-equal for left division");
        }
        matrixDiv(dest, src);
    }

    public static double[][] matrixRightDivGenerate(double[][] dest, double[][] src) throws MadrixException {
        if (dest.length != src.length || dest[0].length != src[0].length) {
            throw new MadrixException("Destination dimensions not-equal for rigth division generation");
        }
        double[][] matrix = dupMatrix(dest);
        matrixDiv(matrix, src);
        return matrix;
    }

    public static void matrixDiv(double[][] dest, double[][] src) throws MadrixException {
        if (dest.length != src.length || dest[0].length != src[0].length) {
            throw new MadrixException("Destination dimensions not-equal for left division");
        }
        for (int i = 0; i < dest.length; i++) {
            MadrixUtils.vectroDiv(dest[i], src[i]);
        }
    }

    public static void vectorSub(double[] dest, double[] subs) throws MadrixException {
        if (dest.length != subs.length) {
            throw new MadrixException("Destination cardinality non-equal to substitutor of vector.");
        }
        for (int i = 0; i < dest.length; i++) {
            dest[i] -= subs[i];
        }
    }

    public static void vectorSum(double[] dest, double[] subs) throws MadrixException {
        if (dest.length != subs.length) {
            throw new MadrixException("Destination cardinality non-equal to the summation of vector.");
        }
        for (int i = 0; i < dest.length; i++) {
            dest[i] += subs[i];
        }
    }

    public static void vectorSum(double[] dest, double value) throws MadrixException {
        for (int i = 0; i < dest.length; i++) {
            dest[i] += value;
        }
    }

    public static void vectroMpl(double[] dest, double[] src) throws MadrixException {
        if (dest.length != src.length) {
            throw new MadrixException("Destination cardinality non-equal to the multiplication of vector.");
        }
        for (int i = 0; i < dest.length; i++) {
            dest[i] *= src[i];
        }
    }

    public static void vectroDiv(double[] dest, double[] src) throws MadrixException {
        if (dest.length != src.length) {
            throw new MadrixException("Destination cardinality non-equal to the division of vector.");
        }
        for (int i = 0; i < dest.length; i++) {
            dest[i] /= src[i];
        }
    }

    public static double vectorPNorm(double[] dest, int p) {
        double norm = 0;
        for (int i = 0; i < dest.length; i++) {
            norm = Math.pow(dest[i], p);
        }
        return Math.pow(norm, 1 / (double) p);
    }

    public static double vectorPDistance(double[] src, double[] des, int p) throws MadrixException {
        if (src.length != des.length) {
            throw new MadrixException("Desctination cardinality non-equal to the src for distance.");
        }
        double norm = 0;
        for (int i = 0; i < des.length; i++) {
            norm += Math.pow(Math.abs(src[i] - des[i]), p);
        }
        return Math.pow(norm, 1 / (double) p);
    }

    public static double vectorDTWDistance(double[] src, double[] des) {
        double[][] dtw = new double[src.length + 1][des.length + 1];
        for (int i = 0; i < src.length; i++) {
            dtw[i][0] = Double.MAX_VALUE;
        }
        for (int j = 0; j < des.length; j++) {
            dtw[0][j] = Double.MAX_VALUE;
        }
        dtw[0][0] = 0;
        for (int i = 1; i < src.length; i++) {
            for (int j = 1; j < des.length; j++) {
                dtw[i][j] = 0;
            }
        }
        for (int i = 0; i < src.length; i++) {
            for (int j = 0; j < des.length; j++) {
                double cost = Math.abs(src[i] - des[j]);
                double mind = Math.min(dtw[i][j + 1], dtw[i + 1][j]);
                mind = Math.min(mind, dtw[i][j]);
                dtw[i + 1][j + 1] = cost + mind;
            }
        }
        return dtw[src.length][des.length];
    }

    public static double vectorDTWDistanceOld(double[] src, double[] des) { //Needs to be checked.
        /*double[] ds = new double[src.length];
        double[] dd = new double[des.length];*/
        double[][] dtw = new double[src.length + 1][des.length + 1];
        for (int i = 0; i < src.length; i++) {
            dtw[i][0] = Double.MAX_VALUE;
        }
        for (int j = 0; j < des.length; j++) {
            dtw[0][j] = Double.MAX_VALUE;
        }
        dtw[0][0] = 0;
        for (int i = 0; i < src.length; i++) {
            for (int j = 0; j < des.length; j++) {
                dtw[i + 1][j + 1] = 0;
            }
        }
        for (int i = 1; i <= src.length; i++) {
            for (int j = 1; j <= des.length; j++) {
                double cost = Math.abs(src[i - 1] - des[j - 1]);
                double mind = Math.min(dtw[i - 1][j], dtw[i][j - 1]);
                mind = Math.min(mind, dtw[i - 1][j - 1]);
                dtw[i][j] = cost + mind;
            }
        }
        return dtw[src.length][des.length];
    }

    public static double matrixDTWPDistance(double[][] src, double[][] des, int p) throws MadrixException {
        if (src[0].length != des[0].length) {
            throw new MadrixException("Matrix DTW pNorm Distance with different dimensions.");
        }
        double[][] dtw = new double[src.length + 1][des.length + 1];
        for (int i = 0; i < src.length; i++) {
            dtw[i][0] = Double.MAX_VALUE;
        }
        for (int j = 0; j < des.length; j++) {
            dtw[0][j] = Double.MAX_VALUE;
        }
        for (int i = 0; i < src.length; i++) {
            for (int j = 0; j < des.length; j++) {
                dtw[i][j] = 0;
            }
        }
        //double[] tmp = new double[src[0].length];
        for (int i = 1; i <= src.length; i++) {
            for (int j = 1; j <= des.length; j++) {
                /*System.arraycopy(src[i - 1], 0, tmp, 0, src[i - 1].length);
                MatrixUtils.vectorSub(tmp, des[i - 1]);
                MatrixUtils.absVector(tmp);
                MatrixUtils.powVector(tmp, p);
                double cost = Math.pow(MatrixUtils.VectorSum(tmp), 1 / (double) p);*/
                double cost = vectorPDistance(src[i - 1], des[j - 1], p);
                double mind = Math.min(dtw[i - 1][j], dtw[i][j - 1]);
                mind = Math.min(mind, dtw[i - 1][j - 1]);
                dtw[i][j] = cost + mind;
            }
        }
        return dtw[src.length][des.length];
    }

    //Matrix permutation.
    public static void permMatrixAtColumn(double[][] dest, int column, int[] perm) throws MadrixException {
        if (perm.length != dest.length) {
            throw new MadrixException("Permutation of Matrix at Column poses error");
        }
        double[] swap = new double[dest.length];
        for (int i = 0; i < dest.length; i++)
            swap[i] = dest[perm[i]][column];
        for (int i = 0; i < dest.length; i++)
            dest[i][column] = swap[i];
    }

    public static void permMatrixColumns(double[][] dest, int[] perm) throws MadrixException {
        if (perm.length != dest[0].length) {
            throw new MadrixException("Permutation of Matrix by Columns poses error");
        }
        double[] swap = new double[dest[0].length];
        for (int i = 0; i < dest.length; i++) {
            for (int j = 0; j < dest[0].length; j++) {
                swap[j] = dest[i][perm[j]];
            }
            System.arraycopy(swap, 0, dest[i], 0, swap.length);
        }
    }

    public static void permMatrixAtRow(double[][] dest, int row, int[] perm) throws MadrixException {
        if (perm.length != dest[0].length) {
            throw new MadrixException("Permutation of Madrix at Row poses error");
        }
        double[] swap = new double[dest[0].length];
        for (int i = 0; i < dest[0].length; i++) {
            swap[i] = dest[row][perm[i]];
        }
        System.arraycopy(swap, 0, dest[row], 0, swap.length);
        /*for (int i = 0; i < dest.length; i++) {
            dest[i] = swap[i];
        }*/
    }

    public static void permMatrixRows(double[][] matrix, int[] perm) throws MadrixException {
        if (perm.length != matrix.length) {
            throw new MadrixException("Permutation of Madrix by Rows poses error.");
        }
        double[][] swap = new double[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            swap[i] = matrix[perm[i]];
        }
        System.arraycopy(swap, 0, matrix, 0, swap.length);
    }

    public static void permVector(double[] vector, int[] perm) {
        double[] swap = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            swap[i] = vector[perm[i]];
        }
        System.arraycopy(swap, 0, vector, 0, vector.length);
        /*for (int i = 0; i < vector.length; i++) {
            vector[i] = swap[i];
        }*/
    }
}
