package org.nlogo.extensions.numanal;

import Jama.Matrix;

public class NAMatrix {
    
    public static Matrix formOuterProduct(Matrix U, Matrix V) {
        // Forms a matrix from the outer product of vectors (1 x n)
        // Matrices U and V.
        return (U.transpose()).times(V);
    }
    
    public static double findMaxAbsElement(Matrix M) {
        // Return the element of M with the largest absolute value by 
        // packing the Matrix into a one dimensional array and then
        // searching for the largest element.
        double[] a = M.getRowPackedCopy();
        double maxElement = 0.0;
        for (double x : a) {
            maxElement = Math.max(maxElement, Math.abs(x));
        }
        return maxElement;
    }
    
    public static Matrix formAbsMatrix(Matrix M) {
        // Return a matrix whose elements are the absolute values of 
        // the elements of M by packing the Matrix into a one dimensional 
        // array, taking the abs of each element, and then converting it 
        // back into a Matrix.
        int nrows = M.getRowDimension();
        double[] a = M.getColumnPackedCopy();
        for (int i = 0; i < a.length; i++) {
            a[i] = Math.abs(a[i]);
        }
        return new Matrix(a, nrows);        
    }
    
        
}