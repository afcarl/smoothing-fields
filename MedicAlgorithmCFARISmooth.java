package edu.jhu.ece.iacl.plugins.dti;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3f;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurface;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.utility.ArrayUtil;

/**
 * 
 * @author Gunnar Atli Sigurdsson, gunnar@jhu.edu
 *
 */
public class MedicAlgorithmCFARISmooth extends ProcessingAlgorithm{
	private ParamVolume basisMixture;
	private ParamVolume basisIndecies;
	private ParamSurface basisVectors;
	private ParamVolume coeffXin;
	private ParamVolume coeffYin;
	private ParamVolume coeffZin;
	private ParamInteger iter;
	
	private ParamVolume smoothDirOut, smoothIndOut, smoothMixOut;
	
	private int rows, cols, slices, components;

	private static final String cvsversion = "$Revision: 1.2 $";
	private static final String revnum = cvsversion.replace("Revision: ", "").replace("$", "") .replace(" ", "");
	private static final String shortDescription = "CFARI Smoothing";
	private static final String longDescription = "Performs N iterations of CFARI Smoothing";


	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(basisIndecies=new ParamVolume("Basis Indecies",VoxelType.INT,-1,-1,-1,-1));
		inputParams.add(basisMixture=new ParamVolume("Mixtures Fraction",VoxelType.FLOAT,-1,-1,-1,-1));
		inputParams.add(basisVectors = new ParamSurface("Basis Directions"));
		inputParams.add(coeffXin=new ParamVolume("CoefficientsX",VoxelType.FLOAT,-1,-1,-1,-1));
		inputParams.add(coeffYin=new ParamVolume("CoefficientsY",VoxelType.FLOAT,-1,-1,-1,-1));
		inputParams.add(coeffZin=new ParamVolume("CoefficientsZ",VoxelType.FLOAT,-1,-1,-1,-1));
		iter = new ParamInteger("Number of Iterations", 1, 1000);
		iter.setValue(4);
		inputParams.add(iter);

		inputParams.setPackage("IACL");
		inputParams.setCategory("DTI.Contrasts");
		inputParams.setLabel("CFARI Smoothing");
		inputParams.setName("CFARI_Smoothing");

		AlgorithmInformation info = getAlgorithmInformation();
		info.setWebsite("");
		info.add(new AlgorithmAuthor("Gunnar Atli Sigurdsson", "gunnar@jhu.edu", "iacl.ece.jhu.edu/gunnar"));
		info.setAffiliation("Johns Hopkins University, Department of Electrical and Computer Engineering");
		info.setDescription(shortDescription);
		info.setLongDescription(shortDescription + longDescription);
		info.setVersion(revnum);
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.ALPHA);
	}


	protected void createOutputParameters(ParamCollection outputParams) {
		smoothDirOut = new ParamVolume("Smoothed Directions",VoxelType.FLOAT,-1,-1,-1,-1);
		smoothDirOut.setName("Smoothed CFARI directions");
		outputParams.add(smoothDirOut);
		
		smoothIndOut = new ParamVolume("Smoothed Direction Approximate Indices",VoxelType.INT,-1,-1,-1,-1);
		smoothIndOut.setName("Smoothed Direction Approximate Indices");
		outputParams.add(smoothIndOut);
		
		smoothMixOut = new ParamVolume("Smoothed Mixtures",VoxelType.FLOAT,-1,-1,-1,-1);
		smoothMixOut.setName("Mixtures after smoothing");
		outputParams.add(smoothMixOut);
	}


	protected void execute(CalculationMonitor monitor) {
		System.out.println(getClass().getCanonicalName()+"\t"+"Getting/Cleaning Data");
		ImageData vecInds = basisIndecies.getImageData();
		String name = vecInds.getName();
		rows = vecInds.getRows();
		cols = vecInds.getCols();
		slices = vecInds.getSlices();
		components = vecInds.getComponents();
		int N = iter.getInt();
		Point3f[] grad;
		int[][][][] dirs;
		float[][][][] mixtures;
		float[][][][] smoothDirFloats;
		Vector3f[][][][] smoothDir;
		Vector3f[][][][] smoothDirPrevious;
		int[][][][] smoothInd;
		float[][][][] smoothMix;
		float[][][][] coeff, coeffX, coeffY, coeffZ;
		
		ImageHeader hdr  = basisIndecies.getImageData().getHeader();
		Point3f resolution = new Point3f();
		resolution.x = hdr.getDimResolutions()[0];
		resolution.y = hdr.getDimResolutions()[1];
		resolution.z = hdr.getDimResolutions()[2];
		System.out.println(getClass().getCanonicalName()+"\t"+"Resolution:" + resolution.x + "," + resolution.y +"," + resolution.z);
		
		
		// this step takes a lot of memory, try/catch
		// TODO is it sufficiently fast to read these coefficients in at every iteration so this can be split into X,Y,Z to save memory?
		try {
			grad = normalizeBasisDirections(basisVectors.getObject());
			basisVectors.dispose();
			dirs = ArrayUtil.getIntArray4d(basisIndecies.getImageData());
			basisIndecies.dispose();
			mixtures = ArrayUtil.getFloatArray4d(basisMixture.getImageData());
			basisMixture.dispose();
			coeffX = ArrayUtil.getFloatArray4d(coeffXin.getImageData());
			coeffXin.dispose();
			coeffY = ArrayUtil.getFloatArray4d(coeffYin.getImageData());
			coeffYin.dispose();
			coeffZ = ArrayUtil.getFloatArray4d(coeffZin.getImageData());
			coeffZin.dispose();
			
			//clean up	
			smoothDir = new Vector3f[rows][cols][slices][components];
			smoothDirPrevious = new Vector3f[rows][cols][slices][components];
			smoothInd = new int[rows][cols][slices][components];
			smoothMix = new float[rows][cols][slices][components];
		}
		catch(Exception exception) {
			System.out.println("Out of memory (or something wrong with inputs)");
			return;
		}
		
		//initialize smoothing
		for(int i = 0; i<rows; i++)
			for(int j = 0; j<cols; j++)
				for(int k = 0; k<slices; k++) {
					
					// get current mixtures and directions
					float[] f = mixtures[i][j][k];
	
					boolean skip = false;
					for(int a=0; a<components; a++) {
						if(Float.isNaN(f[a]))
							skip = true; // faulty voxel
					}
					
					int[] Xfi = dirs[i][j][k];
					
					for(int a=0; a<components; a++) {
						if(Xfi[a] == -1)
							skip = true; // faulty voxel
					}
					
					if (!skip) { // else (0,0,0)
						Vector3f[] Xf = new Vector3f[components];
						for(int a=0; a<components; a++) {
							Xf[a] = new Vector3f(grad[Xfi[a]]);
						}
						
						// do calculations & initialize
						for(int a=0; a<components; a++) {
							smoothDir[i][j][k][a] = Xf[a];
							smoothMix[i][j][k][a] = f[a];
							
							smoothDirPrevious[i][j][k][a] = new Vector3f();
						}
					}
					else {
						for(int a=0; a<components; a++) {
							smoothDir[i][j][k][a] = new Vector3f();
							smoothMix[i][j][k][a] = 0;
							
							smoothDirPrevious[i][j][k][a] = new Vector3f();
						}
						
						
					}
				}
		
		// clean memory
		dirs = null;
		
		
		Point3i p0 = new Point3i();
		Point3i p1 = new Point3i();
		Point3i p2 = new Point3i();
		
		float[] f0, f1, f2;
		Vector3f[] X0, X1, X2; 
		float[] coeff0, coeff1;
		
		float[] f = new float[components];
		Vector3f[] X = new Vector3f[components];
		
		//calculate smoothing
		System.out.println(getClass().getCanonicalName()+"\t"+"Calculating smoothing");
		for(int n = 0; n<N; n++) // iterations
			for(int dim = 0; dim<3; dim++) { // directions
				
				Vector3f[][][][] tmp;
				tmp = smoothDirPrevious;
				smoothDirPrevious = smoothDir;
				smoothDir = tmp;
				
				float[][][][] tmp2;
				tmp2 = mixtures;
				mixtures = smoothMix;
				smoothMix = tmp2;
				
				for(int i = 1; i<rows-1; i++)
					for(int j = 1; j<cols-1; j++)
						for(int k = 1; k<slices-1; k++) {
							p1.set(i, j, k);
							p0.set(p1);
							p2.set(p1);
							
							switch(dim) {
							default:
							case 0:
								p0.x--;
								p2.x++;
								coeff = coeffX;
								break;
							case 1:
								p0.y--;
								p2.y++;
								coeff = coeffY;
								break;
							case 2:
								p0.z--;
								p2.z++;
								coeff = coeffZ;
								break;
							}
							
							// get current mixtures and directions
							f0 = mixtures[p0.x][p0.y][p0.z];
							f1 = mixtures[p1.x][p1.y][p1.z];
							f2 = mixtures[p2.x][p2.y][p2.z];
							
							// remove nans, (set 0)
							for(int a = 0; a < components; a++) {
								if(Float.isNaN(f0[a])) f0[a] = 0;
								if(Float.isNaN(f1[a])) f1[a] = 0;
								if(Float.isNaN(f2[a])) f2[a] = 0;
							}
							
							X0 = smoothDirPrevious[p0.x][p0.y][p0.z];
							X1 = smoothDirPrevious[p1.x][p1.y][p1.z];
							X2 = smoothDirPrevious[p2.x][p2.y][p2.z];
							
							coeff0 = coeff[p0.x][p0.y][p0.z];
							coeff1 = coeff[p1.x][p1.y][p1.z];
							
							// do calculations
							X = oneDimSmooth(f, f0, f1, f2, X0, X1, X2, coeff0, coeff1);
							
							if(Float.isNaN(X[1].x))
								System.out.println("NaN Found.");
			
							for(int a = 0; a < components; a++) {
								smoothDir[i][j][k][a].set(X[a]);
								smoothMix[i][j][k][a] = f[a];
							}
						}
		}
		
		// clean memory
		coeffX = null;
		coeffY = null;
		coeffZ = null;
		coeff = null;
		mixtures = null;
		smoothDirPrevious = null;
		
		
		System.out.println(getClass().getCanonicalName()+"\t"+"CFARI Smooth: Setting up exports.");
		
		// approximate directions
		Vector3f current = new Vector3f();
		for(int i = 1; i<rows-1; i++)
			for(int j = 1; j<cols-1; j++)
				for(int k = 1; k<slices-1; k++)
					for(int a = 0; a < components; a++) {
						int maxi = -1;
						float max = -1;
						for(int b = 0; b < grad.length; b++) {
							
							current.set(grad[b]);
							if(Math.abs(current.dot(smoothDir[i][j][k][a])) > max) {
								max = Math.abs(current.dot(smoothDir[i][j][k][a]));
								maxi = b;
							}
						}
						
						smoothInd[i][j][k][a] = maxi;
				}
		
		
		// free up memory
		grad = null;
		
		// convert from Vector3f to floats for export
		smoothDirFloats = new float[rows][cols][slices][components*3];
		for(int i = 1; i<rows-1; i++)
			for(int j = 1; j<cols-1; j++)
				for(int k = 1; k<slices-1; k++) 
					for(int a = 0; a < components; a++) {
						smoothDirFloats[i][j][k][3*a] = smoothDir[i][j][k][a].x;
						smoothDirFloats[i][j][k][3*a+1] = smoothDir[i][j][k][a].y;
						smoothDirFloats[i][j][k][3*a+2] = smoothDir[i][j][k][a].z;
					}
		
		smoothDir=null; 
		
		// Write outputs

		ImageData out;
		out = (new ImageDataFloat(smoothDirFloats));
		out.setHeader(vecInds.getHeader());
		out.setName(vecInds.getName()+"_smoothDir");
		smoothDirOut.setValue(out);
		smoothDirFloats=null; //dereference to free memory
		
		out = (new ImageDataFloat(smoothMix));
		out.setHeader(vecInds.getHeader());
		out.setName(vecInds.getName()+"_smoothMix");
		smoothMixOut.setValue(out);
		smoothMix=null; //dereference to free memory
		
		out = (new ImageDataInt(smoothInd));
		out.setHeader(vecInds.getHeader());
		out.setName(vecInds.getName()+"_smoothInd");
		smoothIndOut.setValue(out);
		smoothInd=null; //dereference to free memory
	}
	
	
	private Vector3f[] oneDimSmooth(float[] newf1, float[] f0, float[] f1, float[] f2, Vector3f[] X0, Vector3f[] X1, Vector3f[] X2, float[] coeff0, float[] coeff1) {
		
		int N = f1.length; // assume N0 = N1 = N2
		Vector3f[] newX1 = new Vector3f[N];
		for(int i = 0; i < N; i++) {
			newX1[i] = new Vector3f();
		}
		
		float sumF0 = 0;
		float sumF1 = 0;
		float sumF2 = 0;
		for(int i = 0; i < N; i++) {
			sumF0 += f0[i];
			sumF1 += f1[i];
			sumF2 += f2[i];
		}

		
		for(int i = 0; i < N; i++) {
			// newX1[i] = weighted average of all related directions
			
			if( sumF0 + sumF1 + sumF2 < 0.1) { // arbitrary cutoff
				newX1[i].set(X1[i]); // do nothing
				continue;
			}
			
			float[] w0 = new float[N];
			float[] w2 = new float[N];
			
			// find relatives
			int w0i = 0;
			for(int j = i; j < N*N; j = j + N) {
				if(Math.abs(X0[w0i].dot(X1[i])) > 0.7) { //TODO make this parameter
					w0[w0i] = coeff0[j];//sumF0;
				}
				else {
					w0[w0i] = 0.0f;
				}
				w0i++;
			}
			
			int w2i = 0;
			for(int j = i*N; j < i*N+N; j++) {
				if(Math.abs(X2[w2i].dot(X1[i])) > 0.7) {
					w2[w2i] = coeff1[j];//sumF2;
				}
				else {
					w2[w2i] = 0.0f;
				}
				w2i++;
			}
			
			float[] weights = new float[N + 1 + N];
			Vector3f[] dirs = new Vector3f[N + 1 + N];
			
			// concatenate
			//weights[N] = sumF1; // current voxel
			weights[N] = 2*f1[i]; // current voxel
			dirs[N] = X1[i]; // current direction
			for(int j = 0; j < N; j++) {
				dirs[j] = new Vector3f(X0[j]);
				weights[j] = w0[j];
				dirs[N+1+j] = new Vector3f(X2[j]);
				weights[N+1+j] = w2[j];
			}
			
			newX1[i] = weightedDirMean(dirs, weights, X1[i]);
			
			// mixtures
			float contr0 = 0, contr2 = 0;
			for(int j = 0; j < N; j++) {
				contr0 += f0[j]*w0[j];
				contr2 += f2[j]*w2[j];
			}
			
			newf1[i] = contr0/3.0f + f1[i]/3.0f + contr2/3.0f; // FIXME
			newf1[i] = f1[i]; // one step at a time, do this later TODO
		}
			
		return newX1;
	}
	
	private Vector3f weightedDirMean(Vector3f[] dirs, float[] weights, Vector3f orientationDir) {
		// this could be done with clever mathematics, do it simply for speed
		
		Vector3f result = new Vector3f();
		
		// orient dirs to correspond with orientationDir
		for(int i = 0; i < dirs.length; i++) {
			if(dirs[i] == null) {
				return result; // should not happen
			}
			dirs[i].scale(Math.signum(dirs[i].dot(orientationDir))); // flips sign if dot neg
		}
		
		Vector3f tmp = new Vector3f();
		for(int i = 0; i < dirs.length; i++) {
			tmp.set(dirs[i]); 
			tmp.scale(weights[i]);
			result.add(tmp);
		}
		
		if(result.length() < 0.01) {
			return new Vector3f(0, 0, 0);
		}
		else {
			result.normalize();
			return result;
		}
	}
	
	
	
	private Point3f[] normalizeBasisDirections(EmbeddedSurface basisdir){
		Point3f[] dirs = basisdir.getVertexCopy();
		int N = dirs.length;
		float nrm = 0;
		for(int i=0; i<N; i++){
			nrm = 0;
			nrm += dirs[i].x*dirs[i].x;
			nrm += dirs[i].y*dirs[i].y;
			nrm += dirs[i].z*dirs[i].z;
			nrm = (float)Math.sqrt(nrm);
			
//			System.out.println("Before norm dir: " + i + " : " + dirs[i]);
			
			dirs[i].x = (dirs[i].x)/nrm;
			dirs[i].y = (dirs[i].y)/nrm;
			dirs[i].z = (dirs[i].z)/nrm;
			
//			System.out.println("BASIS dir: " + i + " : " + dirs[i]);
		}
		return dirs;
	}
}
