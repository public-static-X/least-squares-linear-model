public class ols_linear_model {
       
    private int N; //number of rows
    private int p; // number of predictors
    private Matrix Beta; // coefficients
    private Matrix standardErrors; //SE of coefficients
    private Matrix T_statistics;

    public ols_linear_model(Matrix X, Matrix Y) {
        if (X.getNumRows() != Y.getNumRows()) throw new IllegalArgumentException("Differing number of rows");
        this.N = X.getNumRows();
        this.p = X.getNumCols();
        
        this.Beta =  ((X.transpose().multiplyMatrix(X)).inverse()).multiplyMatrix(X.transpose()).multiplyMatrix(Y);
        
        Matrix term = Y.subtract(X.multiplyMatrix(Beta));
        double val = term.dotProduct(term);
        double errorVariance = val / (N - p);

        Matrix XtX_inv = ( X.transpose().multiplyMatrix(X) ).inverse();
        Matrix diagMat = Matrix.SQRT( XtX_inv.scalarMultiply(errorVariance) ) ;

        this.standardErrors = diagMat.diag();
        double[][] T_stat = new double[Beta.getNumRows()][2]; 

        for (int i = 0; i < Beta.getNumRows(); i++) {
            T_stat[i][0] = Beta.getElement(i, 0) / standardErrors.getElement(i, 0);
            T_stat[i][1] = 2*T_distribution.cdf(-Math.abs(T_stat[i][0]), N - p - 1);
        }
        T_statistics = new Matrix(T_stat);
    }

    public Matrix coefficients() {
        return Beta;
    }

    public Matrix standardError() {
        return standardErrors;
    }

    @Override
    public String toString() {
        Matrix summary = Beta.columnBind(standardErrors).columnBind(T_statistics);
        StringBuilder names = new StringBuilder();
        names.append( String.format(" %-9s","Beta") );
        names.append( String.format("%-9s","SE") );
        names.append( String.format("%-9s","T stat") );
        names.append( String.format("%-9s","p value") );
        names.append("\n");
        String Summary = summary.toString();
        return names + Summary;
    }

    public static void main(String[] args) {

        Matrix iris = Matrix.read_txt(150,7,"iris_data_sepalLength.txt");
        
        Matrix Y = iris.getColumn(0);
        Matrix X = iris.getColumns(1, 7);
        
        ols_linear_model lm = new ols_linear_model(X, Y);

        System.out.println(lm);
        
    }

}
