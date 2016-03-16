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
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurface;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.utility.ArrayUtil;

import com.cflex.util.lpSolve.LpModel;
import com.cflex.util.lpSolve.LpSolver;

/**
 * 
 * @author Gunnar Atli Sigurdsson, gunnar@jhu.edu
 *
 */
public class MedicAlgorithmCorrespondenceCoefficients extends ProcessingAlgorithm{
	//private ParamSurface cruiseSurf;


	private ParamVolume basisMixture;
	private ParamVolume basisIndecies;
	private ParamSurface basisVectors;
	
	private ParamVolume coeffXout, coeffYout, coeffZout;
	private ParamVolume costXout, costYout, costZout;
	
	private int rows, cols, slices, components;

	private static final String cvsversion = "$Revision: 1.3 $";
	private static final String revnum = cvsversion.replace("Revision: ", "").replace("$", "") .replace(" ", "");
	private static final String shortDescription = "Correspondence Coefficients";
	private static final String longDescription = "Calculates correspondence coefficients.";


	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(basisIndecies=new ParamVolume("Basis Indecies",VoxelType.INT,-1,-1,-1,-1));
		inputParams.add(basisMixture=new ParamVolume("Mixtures Fraction",VoxelType.FLOAT,-1,-1,-1,-1));
		inputParams.add(basisVectors = new ParamSurface("Basis Directions"));

		inputParams.setPackage("IACL");
		inputParams.setCategory("DTI.Contrasts");
		inputParams.setLabel("Correspondence Coefficients");
		inputParams.setName("Correspondence_Coefficients");

		AlgorithmInformation info = getAlgorithmInformation();
		info.setWebsite("");
		info.add(new AlgorithmAuthor("Gunnar Atli Sigurdsson", "gunnar@jhu.edu", "iacl.ece.jhu.edu/gunnar"));
		info.setAffiliation("Johns Hopkins University, Department of Electrical and Computer Engineering");
		info.setDescription(shortDescription);
		info.setLongDescription(shortDescription + longDescription);
		info.setVersion(revnum);
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.BETA);
	}


	protected void createOutputParameters(ParamCollection outputParams) {
		coeffXout = new ParamVolume("CoefficientsX",VoxelType.FLOAT,-1,-1,-1,-1);
		coeffXout.setName("Coefficients in X direction");
		outputParams.add(coeffXout);
		
		coeffYout = new ParamVolume("CoefficientsY",VoxelType.FLOAT,-1,-1,-1,-1);
		coeffYout.setName("Coefficients in Y direction");
		outputParams.add(coeffYout);
		
		coeffZout = new ParamVolume("CoefficientsZ",VoxelType.FLOAT,-1,-1,-1,-1);
		coeffZout.setName("Coefficients in Z direction");
		outputParams.add(coeffZout);
		
		costXout = new ParamVolume("CostX",VoxelType.FLOAT,-1,-1,-1,1);
		costXout.setName("Cost in X direction");
		outputParams.add(costXout);

		costYout = new ParamVolume("CostY",VoxelType.FLOAT,-1,-1,-1,1);
		costYout.setName("Cost in Y direction");
		outputParams.add(costYout);
		
		costZout = new ParamVolume("CostZ",VoxelType.FLOAT,-1,-1,-1,1);
		costZout.setName("Cost in Z direction");
		outputParams.add(costZout);
	}


	protected void execute(CalculationMonitor monitor) {
		System.out.println(getClass().getCanonicalName()+"\t"+"Calculating coefficients");
		ImageData vecInds = basisIndecies.getImageData();
		String name = vecInds.getName();
		rows = vecInds.getRows();
		cols = vecInds.getCols();
		slices = vecInds.getSlices();
		components = vecInds.getComponents();
		Point3f[] grad;
		int[][][][] dirs;
		float[][][][] mixtures;
		float[][][] costX;
		float[][][] costY;
		float[][][] costZ;
		float[][][][] coeffX;
		float[][][][] coeffY;
		float[][][][] coeffZ;
		
		ImageHeader hdr  = basisIndecies.getImageData().getHeader();
		Point3f resolution = new Point3f();
		resolution.x = hdr.getDimResolutions()[0];
		resolution.y = hdr.getDimResolutions()[1];
		resolution.z = hdr.getDimResolutions()[2];
		System.out.println(getClass().getCanonicalName()+"\t"+"Resolution:" + resolution.x + "," + resolution.y +"," + resolution.z);
		
		
		// this step takes alot of memory, try/catch
		// TODO calculate X,Y,Z separately to save memory
		try {
			grad = normalizeBasisDirections(basisVectors.getObject());
			dirs = ArrayUtil.getIntArray4d(basisIndecies.getImageData());
			mixtures = ArrayUtil.getFloatArray4d(basisMixture.getImageData());
			
			//clean up
			basisIndecies.dispose();
			basisMixture.dispose();
			basisVectors.dispose();
	
			costX = new float[rows][cols][slices];
			costY = new float[rows][cols][slices];
			costZ = new float[rows][cols][slices];
			
			coeffX = new float[rows][cols][slices][components*components];
			coeffY = new float[rows][cols][slices][components*components];
			coeffZ = new float[rows][cols][slices][components*components];
		}
		catch(Exception exception) {
			System.out.println("Out of memory (or something wrong with inputs)");
			return;
		}
		
		Point3i p = new Point3i();
		Point3i p2 = new Point3i();
		Vector3f[] Xf = new Vector3f[components];
		Vector3f[] Xg = new Vector3f[components];
		
		for(int i=0; i < components; i++) {
			Xf[i] = new Vector3f();
			Xg[i] = new Vector3f();
		}
		
		//calculate coefficients
		for(int i = 0; i<rows-1; i++)
			for(int j = 0; j<cols-1; j++)
				for(int k = 0; k<slices-1; k++)
					for(int dim = 0; dim<2; dim++) {
						// for each dimension
						p.set(i, j, k);
						p2.set(p);
						
						switch(dim) {
						case 0:
							p2.x++;
							break;
						case 1:
							p2.y++;
							break;
						case 2:
							p2.z++;
							break;
						}
						
						// get current mixtures and directions
						float[] f = mixtures[p.x][p.y][p.z];
						float[] g = mixtures[p2.x][p2.y][p2.z];
		
						boolean skip = false;
						for(int a=0; a<components; a++) {
							if(Float.isNaN(f[a]) || Float.isNaN(g[a]))
								skip = true; // faulty voxel
						}
						
						int[] Xfi = dirs[p.x][p.y][p.z];
						int[] Xgi = dirs[p2.x][p2.y][p2.z];
						
						for(int a=0; a<components; a++) {
							if(Xfi[a] == -1 || Xgi[a] == -1)
								skip = true; // faulty voxel
						}
						
						if (!skip) {
							for(int a=0; a<components; a++) {
								Xf[a].set(grad[Xfi[a]]);
								Xg[a].set(grad[Xgi[a]]);
							}
							
							
							// finally, calculate coefficients and cost
							switch(dim) {
							case 0:
								costX[p.x][p.y][p.z] = oneDimCoeffs(f, g, Xf, Xg, coeffX[p.x][p.y][p.z]);
								break;
							case 1:
								costY[p.x][p.y][p.z] = oneDimCoeffs(f, g, Xf, Xg, coeffY[p.x][p.y][p.z]);
								break;
							case 2:
								costZ[p.x][p.y][p.z] = oneDimCoeffs(f, g, Xf, Xg, coeffZ[p.x][p.y][p.z]);
								break;
							}
							
						}
					}
		
		// free up memory memory
		grad = null;
		mixtures = null;
		dirs = null;
		
		// Write outputs
		System.out.println(getClass().getCanonicalName()+"\t"+"Correspondence Coefficients: Setting up exports.");

		ImageData out;
		out = (new ImageDataFloat(costX));
		out.setHeader(vecInds.getHeader());
		out.setName(vecInds.getName()+"_costX");
		costXout.setValue(out);
		costX=null; //dereference to free memory

		out = (new ImageDataFloat(costY));
		out.setHeader(vecInds.getHeader());
		out.setName(vecInds.getName()+"_costY");
		costYout.setValue(out);
		costY=null; //dereference to free memory
		
		out = (new ImageDataFloat(costZ));
		out.setHeader(vecInds.getHeader());
		out.setName(vecInds.getName()+"_costZ");
		costZout.setValue(out);
		costZ=null; //dereference to free memory
		
		out = (new ImageDataFloat(coeffX));
		out.setHeader(vecInds.getHeader());
		out.setName(vecInds.getName()+"_coeffX");
		coeffXout.setValue(out);
		coeffX=null; //dereference to free memory
		
		out = (new ImageDataFloat(coeffY));
		out.setHeader(vecInds.getHeader());
		out.setName(vecInds.getName()+"_coeffY");
		coeffYout.setValue(out);
		coeffY=null; //dereference to free memory
		
		out = (new ImageDataFloat(coeffZ));
		out.setHeader(vecInds.getHeader());
		out.setName(vecInds.getName()+"_coeffZ");
		coeffZout.setValue(out);
		coeffZ=null; //dereference to free memory
	}
	
	
	
	/**
	 * 
	 * Calculates coefficients between two voxels
	 * 
	 * @param f mixtures of first voxel
	 * @param g mixtures of second voxel
	 * @param Xf directions of first voxel
	 * @param Xg directions of second voxel
	 * @param coefficients the output coefficients 
	 * @return The cost of the transformation
	 */
	private float oneDimCoeffs(float[] f, float[] g, Vector3f[] Xf, Vector3f[] Xg, float[] coeff) {
		int Nf = f.length;
		int Ng = g.length;
		
		// Normalize mixtures;
		float[] fnorm = new float[Nf];
		float[] gnorm = new float[Ng];
		
		float sumF = 0;
		float sumG = 0;
		
		for(int i = 0; i < Nf; i++) {
			sumF += f[i];
		}
		
		for(int i = 0; i < Ng; i++) {
			sumG += g[i];
		}
		
		for(int i = 0; i < Nf; i++) {
			fnorm[i] = f[i]/sumF;
		}
			
		for(int i = 0; i < Ng; i++) {
			gnorm[i] = g[i]/sumG;
		}
		
		
		// create f matrix for linear program
		double[][] F = new double[Ng][Ng*Nf+1]; // +1 for this silly linprog solver
		for(int i = 0; i < Ng*Nf; i++) {
			int elem = i/Ng;
			int row = i - Ng*(elem);
			
			F[row][i+1] = 1; 
		}
		
		
		// calculate movement penalties
		int iD = 1;
		double chord;
		double arclength;
		Vector3f vec1 = new Vector3f();
		Vector3f vec2 = new Vector3f();
		double[] D = new double[Nf*Ng+1];
		
		for(int i = 0; i<Nf; i++) // index for first voxel
			for(int j = 0; j<Ng; j++) { // index for second voxel
				vec1.sub(Xf[i], Xg[j]);
				Xf[i].negate();
				vec2.sub(Xf[i], Xg[j]);
				
				if(vec2.lengthSquared() < vec1.lengthSquared()) {
					vec1.set(vec2);
				}
				
				chord = vec1.length();
				arclength = 2*Math.asin(chord/2);
				D[iD] = arclength*arclength;
				
				iD++;
			}
		
		double[][] S = new double[Ng][Ng*Nf+1];
		
		int iS = 1;
		for(int i = 0; i < Ng; i++) {
			// for all rows
			for(int j=0; j<Nf; j++) {
				S[i][iS++] = 1;
			}
		}
		
		// Solve linear problem using F, S and D to get the coefficients
        try {
            LpModel lpIn = new LpModel(0, Ng*Nf);

            lpIn.setObjFn(D); // Objective function
            for (int i = 0; i<Ng; i++) {
                lpIn.addConstraint(F[i], LpModel.EQ, gnorm[i]); // equality constraint
                lpIn.setLowerBound(i+1, 0); // lower bound
                lpIn.addConstraint(S[i], LpModel.LE, fnorm[i]); // upper bound, better performance
                //lpIn.addConstraint(S[i], LpModel.EQ, fnorm[i]);
            }

            LpSolver lpSolve = new LpSolver(lpIn);
            int flag = lpSolve.solve();
            
            if (flag!=0)
                return 0.0f; // error
            
            for (int i = 0; i<Nf*Ng; i++)
                coeff[i] = (float) lpIn.getResult(i+1);
            
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        
        // find cost
        float cost = 0.0f;
        
        for(int i=0; i < Nf*Ng; i++) {
        	cost += coeff[i]*D[i];        	
        }
        
//        // convert to shared mixture values
//        for(int i=0; i < Nf*Ng; i++) {
//        	coeff[i] = coeff[i]*fnorm[i/Ng];
//        }
       		
		return cost;
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
