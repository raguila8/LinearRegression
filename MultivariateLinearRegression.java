import java.io.*;
import java.util.*;
import javax.swing.*;
import org.math.plot.*;
import static org.math.io.files.ASCIIFile.*;
import static org.math.io.parser.ArrayString.*;

/*************************************************************************
* 
*
*
*
*
*
*
*
*
****************************************************************************/


public class MultivariateLinearRegression {
    
    //Vector<Vector<Double>> trainingSet;
		private double[][] trainingSet;
    //Vector<Double> theta = new Vector<Double>();
		private double[] theta;
    private double alpha = 0.000001;
    private double tol = 1e-12;
    private int maxiter = 10000000;
		private int dispiter = 1000;
		private int iter = 0;
		// Number of training examples
		private int trainingExamples;
    // Number of features
    private int features;
		private double[] cost;
		private double[] iterations;

		//Plot3DPanel plot;
    

		// trainingSet is a 2D array where each row corresponds to a particular 
		// feature and its columns are its values. The last row correspnds to the
		// output.
    public MultivariateLinearRegression(double[][] trainingSet) {
	    this.trainingSet = trainingSet;
			this.features = this.trainingSet.length - 1;
			this.trainingExamples = this.trainingSet[0].length;
		}
	
		public void plotData() {
			if (this.features >= 3) {
				System.err.println("Error: Cannot print plot with 3 or more independent variables");
				return;
			} else {
				if (this.features == 1) {
					Plot2DPanel plot = new Plot2DPanel();
					this.plot.addScatterPlot("X-Y", this.trainingSet[0], this.trainingSet[1]);
					JFrame frame = new JFrame("Final X-Y Data");
					frame.setContentPane(this.plot);
					frame.setSize(600, 600);
					frame.setVisible(true);
				} else if (this.features == 2) {
					Plot3DPanel plot = new Plot3DPanel();
					this.plot.addScatterPlot("X-Y-Z", this.trainingSet[0], this.trainingSet[1], this.trainingSet[2]);
					JFrame frame = new JFrame("Final X-Y-Z Data");
					frame.setContentPane(this.plot);
					frame.setSize(600, 600);
					frame.setVisible(true);
				}
			}
		}

		public void addTrendLine() {
			double[] yEnd = new double[this.trainingExamples];
			for (int i = 0; i < this.trainingExamples; i++) {
				//yEnd[i] = 
			}
		}
    
    public double getTheta(int index) {
        //return this.theta.get(index);
				return this.theta[index];
    }

		public int getMaxiter() {
			return this.maxiter;
		} 

		public int getDispiter() {
			return this.dispiter;
		}

		public double getAlpha() {
			return this.alpha;
		}

		public double getTol() {
			return this.tol;
		}

		public void setMaxiter(int maxiter) {
			this.maxiter = maxiter;
		}

		public void setDispiter(int dispiter) {
			this.dispiter = dispiter;
		}

		public void setAlpha(double alpha) {
			this.alpha = alpha;
		} 

		public void setTol(double tol) {
			this.tol = tol;
		}

    // The input is of size features + 1. The zeroth "feature" is 1.
    public double hypothesisFunction(double[] input) {
        double sum = 0;
        for (int i = 0; i < features + 1; i++) {
            //sum += input.get(i) * this.theta.get(i);
						sum += input[i] * this.theta[i];
        }
        return sum;
    }

		private double costFunction() {
			double sum = 0;

			double[] input = new double[this.features + 1];
			input[0] = 1.0;
			for (int j = 0; j < this.trainingExamples; j++) {
				for (int i = 0; i < this.features; i++) {
					input[i + 1] = this.trainingSet[i][j];
					if (i == this.features - 1) {
							sum += Math.pow((this.hypothesisFunction(input) - this.trainingSet[this.features][j]), 2);
					}
				}					
			}
			return ((double) 1/(2 * this.trainingExamples)) * sum;
		}
    
    
    private double deriveTheta(int index) {
        double sum = 0;

        //Vector<Double> input = new Vector<Double>();
        //input.add(1.0);
			double[] input = new double[this.features + 1];
			for (int j = 0; j < this.trainingExamples; j++) {
				for (int i = 0; i < this.features; i++) {
					input[i + 1] = this.trainingSet[i][j];
					if (i == this.features - 1) {
						if (index != 0) {
							sum += (this.hypothesisFunction(input) - this.trainingSet[this.features][j]) 
									* this.trainingSet[index - 1][j];
						} else {
							sum += (this.hypothesisFunction(input) - this.trainingSet[this.features][j]); 
						}
					}
				}					
			}
			return sum;
    }
    
    public void execute() {
        int count = 0;
				this.theta = new double[this.features + 1];
        for (int i = 0; i < this.features + 1; i++) {
						this.theta[i] = 0.0;
        }
				this.cost = new double[(int) Math.ceil((double)this.maxiter / this.dispiter) + 1];
				this.iterations = new double[this.cost.length];

        do {
            //double[] temp = new double[this.theta.size()];
						double[] temp = new double[this.theta.length];

            for (int i = 0; i < temp.length; i++) {
                //temp[i] = this.theta.get(i);
								temp[i] = this.theta[i];

                //temp[i] -= this.alpha * ((double) 1/this.trainingSet.get(i).size()) * this.deriveTheta(i);
								temp[i] -= this.alpha * ((double) 1/this.trainingSet[i].length) * this.deriveTheta(i);

            }
            for (int i = 0; i < temp.length; i++) {
                //this.theta.set(i, temp[i]);
								this.theta[i] = temp[i];
            }

						if (this.iter % this.dispiter == 0) {
							this.cost[count] = this.costFunction();
							this.iterations[count] = this.iter;
							if (count >= 2) {
								if ((this.cost[count] - this.cost[count - 1]) > 0) {
            	    this.alpha = (double) this.alpha / (1 + ((double) (this.iter - (this.iter / 1000))/ 10000));
                }
							//System.out.println(Math.abs(this.cost[count] - this.cost[count - 1]) < this.tol);
							if (Math.abs(this.cost[count] - this.cost[count - 1]) < this.tol) {
									this.iter = count;
									break;
							}
						}

							count++;
						}		
            
            this.iter++;
        } while (this.iter <= this.maxiter);
    }

		public void printConvergence() {
			// plot the convergednce data
			Plot2DPanel convPlot = new Plot2DPanel();

			if (this.iter < this.maxiter) {
				double[] temp0 = new double[this.iter];
				double[] temp1 = new double[this.iter];
				for (int i = 0; i < this.iter; i++) {
					temp0[i] = this.cost[i];
					temp1[i] = this.iterations[i];
				}
				this.cost = temp0;
				this.iterations = temp1;
			}
			// add a line to the PlotPanel
			convPlot.addLinePlot("Cost Function", this.iterations, this.cost);
			
			// put the PlotPanel in a JFrame, as a JPanel
			JFrame frame = new JFrame("Convergence of the Cost Function");
			//frame.setDefaultCloseOperation(EXIT_ON_CLOSE);
			frame.setContentPane(convPlot);
			frame.setSize(600, 600);
			frame.setVisible(true);
		}
    
    public static void main(String[] args) {
        Scanner in = new Scanner(System.in);
        int features = in.nextInt();
        int trainingExamples = in.nextInt();
				double[][] trainingSet = new double[features + 1][trainingExamples];
        //Vector<Vector<Double>> trainingSet = new Vector<Vector<Double>>();
        //for (int i = 0; i < features + 1; i++) {
            //trainingSet.add(new Vector<Double>());
        //}
        for (int i = 0; i < trainingExamples; i++) {
            for (int j = 0; j < features + 1; j++) {
								trainingSet[j][i] = in.nextDouble();
                //trainingSet.get(j).add(in.nextDouble());
            }
        }
        int numberOfInputs = in.nextInt();
        //Vector<Vector<Double>> input = new Vector<Vector<Double>>();
        double[][] input = new double[numberOfInputs][features + 1];
        for (int i = 0; i < numberOfInputs; i++) {
            //input.add(new Vector<Double>());
            input[i][0] = 1.0;
            for (int j = 1; j < features + 1; j++) {
                //input.get(i).add(in.nextDouble());
								input[i][j] = in.nextDouble();
            }
        }
        
        MultivariateLinearRegression s = new MultivariateLinearRegression(trainingSet);
				//s.plotData();
        s.execute();
        double[] output = new double[numberOfInputs];
        for (int i = 0; i < numberOfInputs; i++) {
            //output[i] = s.hypothesisFunction(input.get(i));
						output[i] = s.hypothesisFunction(input[i]);
            System.out.println(output[i]);
        }
				//s.addTrendLine();
				s.printConvergence();
    }
}
