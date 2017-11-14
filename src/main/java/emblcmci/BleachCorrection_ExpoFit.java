package emblcmci;

/** Bleach Correction by Fitting Exponential Decay function.
 *  Kota Miura (miura@embl.de)
 *
 * 2D and 3D time series corrected by fitting exponential decay function.
 *
 * Copyright © 2010 Kota Miura
 * License: GPL 2
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 2
 * as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.gui.Plot;
import ij.measure.CurveFitter;
import ij.plugin.frame.Fitter;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.gui.GenericDialog;
import ij.util.Tools;
import java.awt.Color;

public class BleachCorrection_ExpoFit {
	ImagePlus imp;
	boolean is3DT = false;
	Roi curROI = null;
	int fitStartFrame;
	int fitEndFrame;
	Plot fitplot;
	/**
	 * @param imp
	 */
	public BleachCorrection_ExpoFit(ImagePlus imp) {
		super();
		this.imp = imp;
	}

	public BleachCorrection_ExpoFit(ImagePlus imp, Roi curROI) {
		super();
		this.imp = imp;
		this.curROI = curROI;
	}

	/** Fit the mean intensity time series of given ImagePlus in this class.
	 * fitStartFrame: first fame number for fitting, starting from 1. 
	 * fitEndFrame: last framenumber to be fitted
	 * fit equation is 11, parameter from
	 * 		http://rsb.info.nih.gov/ij/developer/api/constant-values.html#ij.measure.CurveFitter.STRAIGHT_LINE
	 * @return an instance of CurveFitter
	 */
	public CurveFitter dcayFitting(){
		ImageProcessor curip;
		ImageStatistics imgstat;
		int fitframes = fitEndFrame - fitStartFrame + 1;
		double[] xA = new double[imp.getStackSize()];
		double[] yA = new double[imp.getStackSize()];
		double[] fitxA = new double[fitframes];
		double[] fityA = new double[fitframes];

		if (curROI == null) curROI = new Roi(0, 0, imp.getWidth(), imp.getHeight());

		int counter = 0;
		for (int i = 0; i < imp.getStackSize(); i++){
			curip = imp.getImageStack().getProcessor(i + 1);
			curip.setRoi(curROI);
			imgstat = curip.getStatistics();
			xA[i] = i;
			yA[i] = imgstat.mean;
			if ((i >= fitStartFrame - 1) && (i <= fitEndFrame - 1)){
				fitxA[counter] = i;
				fityA[counter] = imgstat.mean;
				counter++;
			}
		}
		CurveFitter cf = new CurveFitter(fitxA, fityA);
		double firstframeint = fityA[0];
		double lastframeint = fityA[fityA.length - 1];
		double guess_a = firstframeint - lastframeint;
		if (guess_a <= 0){
			IJ.error("This sequence seems to be not decaying");
			return null;
		}
		double guess_c = lastframeint;
		double maxiteration = 2000;
		double NumRestarts = 2;
		double errotTol = 10;
		double[] fitparam = {-1*guess_a, -0.0001, guess_c, maxiteration, NumRestarts, errotTol};

		cf.setInitialParameters(fitparam);
		cf.doFit(CurveFitter.EXP_WITH_OFFSET); 
		//Fitter.plot(cf);
		IJ.log(cf.getResultString());
		fitplot = plotter(cf, xA, yA);
		return cf;
	}

	/** Curve fitting for 3D time series is done with average intensity value for
	 * each time point (per time point intensity mean is used, so the fitted points = time point, not slice number)
	 *
	 * @param zframes
	 * @param tframes
	 * @return
	 */
	public CurveFitter decayFitting3D(int zframes, int tframes){
		ImageProcessor curip;
		ImageStatistics imgstat;
		double[] xA = new double[tframes];
		double[] yA = new double[tframes];
		double curStackMean = 0.0;
		if (curROI == null) curROI = new Roi(0, 0, imp.getWidth(), imp.getHeight());
		for (int i = 0; i < tframes; i++){
			curStackMean = 0.0;
			for (int j = 0; j < zframes; j++){
				curip = imp.getImageStack().getProcessor(i * zframes + j + 1);
				curip.setRoi(curROI);
				imgstat = curip.getStatistics();
				curStackMean += imgstat.mean;
			}
			curStackMean /= zframes;
			xA[i] = i;
			yA[i] =	curStackMean;
		}
		CurveFitter cf = new CurveFitter(xA, yA);
		double firstframeint = yA[0];
		double lastframeint = yA[yA.length - 1];
		double guess_a = firstframeint - lastframeint;
		if (guess_a <= 0){
			IJ.error("This sequence seems to be not decaying");
			return null;
		}
		double guess_c = lastframeint;
		double maxiteration = 2000;
		double NumRestarts = 2;
		double errotTol = 10;
		double[] fitparam = {-1 * guess_a, -0.0001, guess_c, maxiteration, NumRestarts, errotTol};

		cf.setInitialParameters(fitparam);

		cf.doFit(11); 
		Fitter.plot(cf);
		IJ.log(cf.getResultString());
		return cf;
	}


	/** calculate estimated value from fitted "Exponential with Offset" equation
	 *
	 * @param a  magnitude (difference between max and min of curve)
	 * @param b  exponent, defines degree of decay
	 * @param c  offset.
	 * @param x  timepoints (or time frame number)
	 * @return estimate of intensity at x
	 */
	public double calcExponentialOffset(double a, double b, double c, double x){
		return (a * Math.exp(-b*x) + c);
	}

	void setFitRange(){
		GenericDialog gd = new GenericDialog("Set Fitting Start and End frame numbers");
		gd.addNumericField("Start fame: ", fitStartFrame, 0);
		gd.addNumericField("End frame: ", fitEndFrame, 0);
		gd.showDialog();
		if (gd.wasCanceled()) return;
		if (fitEndFrame <= fitStartFrame){
			IJ.showMessage("End frame should be larger than the start frame!");
			return;
		}
		int newfitStartFrame = (int)gd.getNextNumber();
		int newfitEndFrame = (int)gd.getNextNumber();
		if (newfitStartFrame < 1)
			fitStartFrame = 1;
		else
			fitStartFrame = newfitStartFrame;
		if ((newfitEndFrame <= fitEndFrame) && (newfitEndFrame > 1))
			fitEndFrame = newfitEndFrame;
	}

	/** Plot results. 
	 * modified plot method in Fitter class
	 * see https://github.com/imagej/imagej1/blob/master/ij/plugin/frame/Fitter.java
	 */
	Plot plotter(CurveFitter cf, double[] xA, double[] yA){
		double[] x = cf.getXPoints();
		double[] y = cf.getYPoints();		
		if (cf.getParams().length<cf.getNumParams()) {
			Plot plot = new Plot(cf.getFormula(),"X","Y",x,y);
			plot.setColor(Color.BLUE);
			plot.addLabel(0.02, 0.1, cf.getName());
			plot.addLabel(0.02, 0.2, cf.getStatusString());
			plot.show();
			return null;
		}
		int npoints = 100;
		if (npoints<x.length)
			npoints = x.length; //or 2*x.length-1; for 2 values per data point
		if (npoints>1000)
			npoints = 1000;
		double[] a = Tools.getMinMax(x);
		double xmin=a[0], xmax=a[1];
		//if (eightBitCalibrationPlot) {
			//npoints = 256;
			//xmin = 0;
			//xmax = 255;
		//}
		a = Tools.getMinMax(y);
		double ymin=a[0], ymax=a[1]; //y range of data points
		float[] px = new float[npoints];
		float[] py = new float[npoints];
		double inc = (xmax-xmin)/(npoints-1);
		double tmp = xmin;
		for (int i=0; i<npoints; i++) {
			px[i]=(float)tmp;
			tmp += inc;
		}
		double[] params = cf.getParams();
		for (int i=0; i<npoints; i++)
			py[i] = (float)cf.f(params, px[i]);
		a = Tools.getMinMax(py);
		double[] dataXMinMax = Tools.getMinMax(xA);
		double[] dataYMinMax = Tools.getMinMax(yA);
		double dataRange = ymax - ymin;
		//ymin = Math.max(ymin - dataRange, Math.min(ymin, a[0])); //expand y range for curve, but not too much
		//ymax = Math.min(ymax + dataRange, Math.max(ymax, a[1]));
		double plotxMin = -1;//dataXMinMax[0] * 0.9;
		double plotxMax = dataXMinMax[1] * 1.1;
		double plotyMin = dataYMinMax[0] * 0.9;
		double plotyMax = dataYMinMax[1] * 1.1;
		//Plot plot = new Plot(cf.getFormula(),"X","Y",px,py);
		Plot plot = new Plot(cf.getFormula(),"Frame","Mean Intensity");
		plot.setLimits(plotxMin, plotxMax, plotyMin, plotyMax);
		plot.setColor(Color.BLUE);
		plot.setLineWidth(2);
		plot.addPoints(px,py,Plot.LINE);
		//plot.setLimits(xmin, xmax, ymin, ymax);
		plot.setColor(Color.RED);
		plot.setLineWidth(1);
		plot.addPoints(xA, yA, Plot.CIRCLE);
		plot.setColor(Color.BLUE);

		StringBuffer legend = new StringBuffer(100);
		legend.append(cf.getName()); legend.append('\n');
		legend.append(cf.getFormula()); legend.append('\n');
        double[] p = cf.getParams();
        int n = cf.getNumParams();
        char pChar = 'a';
        for (int i = 0; i < n; i++) {
			legend.append(pChar+" = "+IJ.d2s(p[i],5,9)+'\n');
			pChar++;
        }
		legend.append("R^2 = "+IJ.d2s(cf.getRSquared(),4)); legend.append('\n');
		plot.addLabel(0.02, 0.1, legend.toString());
		plot.setColor(Color.BLUE);
		return plot;//lot.show();				
	}

	/** does both decay fitting and bleach correction.
	 *
	 */
	public void core(){
		int[] impdimA = imp.getDimensions();
		IJ.log("slices"+Integer.toString(impdimA[3])+"  -- frames"+Integer.toString(impdimA[4]));
		//IJ.log(Integer.toString(imp.getNChannels())+":"+Integer.toString(imp.getNSlices())+":"+ Integer.toString(imp.getNFrames()));
		int zframes = impdimA[3];
		int tframes = impdimA[4];
		if (zframes > 1 && tframes > 1){	// if slices and frames are both more than 1
			is3DT =true;
			if ((zframes * tframes) != imp.getStackSize()){
				IJ.showMessage("slice and time frames do not match with the length of the stack. Please correct!");
				return;
			}
		} else {
			fitStartFrame = 1;
			if (zframes > 1) 
				fitEndFrame = zframes;
			else 
				fitEndFrame = tframes;
		}
		CurveFitter cf;
		if (is3DT) cf = decayFitting3D(zframes, tframes);
		else { 
			setFitRange();
			cf = dcayFitting();
		}
		double[] respara = cf.getParams();
		double res_a = respara[0];
		double res_b = respara[1];
		double res_c = respara[2];
		double ratio = 0.0;
		ImageProcessor curip;
		System.out.println(res_a + "," + res_b + "," + res_c);
		double[] outxA, outyA;
		if (is3DT){
			for (int i = 0; i < tframes; i++){
				for (int j = 0; j < zframes; j++){
					curip = imp.getImageStack().getProcessor(i * zframes + j + 1);
					ratio = calcExponentialOffset(res_a, res_b, res_c, 0.0) / calcExponentialOffset(res_a, res_b, res_c, (double) (i + 1));
					curip.multiply(ratio);
				}
			}
		} else {
			outxA = new double[imp.getStackSize()];
			outyA = new double[imp.getStackSize()];
			for (int i = 0; i < imp.getStackSize(); i++){
				curip = imp.getImageStack().getProcessor(i+1);
				ratio = calcExponentialOffset(res_a, res_b, res_c, 0.0) / calcExponentialOffset(res_a, res_b, res_c, (double) (i + 1));
				curip.multiply(ratio);
				outxA[i] = i;
				outyA[i] = curip.getStatistics().mean;
			}
			fitplot.setColor(Color.GREEN);
			fitplot.addPoints(outxA, outyA, Plot.CIRCLE);
			fitplot.show();
		}
	}



}
