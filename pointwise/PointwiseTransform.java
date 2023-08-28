public class PointwiseTransform extends Object {

	/**
	* Question 2.1 Contrast reversal
	*/
	static public ImageAccess inverse(ImageAccess input) {
		int nx = input.getWidth();
		int ny = input.getHeight();
		ImageAccess output = new ImageAccess(nx, ny);
		double value = 0.0;
		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++) {
			value = input.getPixel(x, y);
			value = 255 - value;
			output.putPixel(x, y, value);
		}
		return output;	
	}

	/**
	* Question 2.2 Stretch normalized constrast
	*/
	static public ImageAccess rescale(ImageAccess input) {
	    // Get Image with / height
		int nx = input.getWidth();
		int ny = input.getHeight();

		// Calculate alpha coefficient
		double max = input.getMaximum();
		double min = input.getMinimum();
	    double alpha = 255.0 / (max - min);
	    
		// Calculate beta coefficient
	    double beta = min;
	
	    ImageAccess output = new ImageAccess(nx, ny);
	
	    double value = 0.0;
	    for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++) {
			value = alpha * (input.getPixel(x, y) - beta);
			output.putPixel(x, y, value);
		}
		return output;	
	}

	/**
	* Question 2.3 Saturate an image
	*/
	static public ImageAccess saturate(ImageAccess input) {
		// Get width / Height
		int nx = input.getWidth();
		int ny = input.getHeight();
		
		ImageAccess output = new ImageAccess(nx, ny);
	
		// If pixel value is more than 10000, saturate the pixel  
		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++) {
			if (input.getPixel(x, y) > 10000) {
				output.putPixel(x, y, 10000);
			} else {
				output.putPixel(x, y, input.getPixel(x, y));
			}
		}
		return output;
	}
	
	/**
	* Question 4.1 Maximum Intensity Projection
	*/
	static public ImageAccess zprojectMaximum(ImageAccess[] zstack) {
		// Get width / height and Z
    	int nx = zstack[0].getWidth();
		int ny = zstack[0].getHeight();
		int nz = zstack.length;
		
    	ImageAccess output = new ImageAccess(nx, ny);
	  
    	double value_max = 0.0;
    	for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++) {
    	  	value_max = zstack[0].getPixel(x, y);
    	  	// Iterate between Z axis
    	  	for (int z=0; z<nz; z++){
				if (value_max < zstack[z].getPixel(x, y)) {
    	      	value_max = zstack[z].getPixel(x, y);
				}
    	  	}

			output.putPixel(x, y, value_max);
    	}
		return output;	
	}

	/**
	* Question 4.2 Z-stack mean
	*/
	static public ImageAccess zprojectMean(ImageAccess[] zstack) {
		int nx = zstack[0].getWidth();
		int ny = zstack[0].getHeight();
		int nz = zstack.length;
		
		ImageAccess output = new ImageAccess(nx, ny);
			
    	double z_accumulator = 0.0;
    	for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++) {
			z_accumulator = 0.0;
    		for (int z=0; z<nz; z++){
    	    	z_accumulator += zstack[z].getPixel(x, y); 
		  	}
			output.putPixel(x, y, z_accumulator/nz);
    	}
	
		return output;	
	}

}
