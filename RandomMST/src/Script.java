import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class Script {

	public Script() {}
	
	public static void main(String[] args) throws IOException {
		int[] dims = {1,2,3,4};
		int[] points = {16384, 32768, 65536, 131072}; //128, 256, 512, 1024, 2048, 4096, 8192, 
		int flag = 1;
		int n_trials = 1;
		
//		File file = new File("../results.csv");
//		if (!file.exists()) 
//			file.createNewFile();
//		BufferedWriter bw = new BufferedWriter(new FileWriter(file));
		
		double avg_weight = -1;		
		for (int j = 0; j < points.length; j++) {
			for (int i = 0; i < dims.length; i++) {
				avg_weight = RandMST.runExperiment(flag, points[j], n_trials, dims[i]);
				
//				bw.write(String.valueOf(avg_weight));
//				bw.write(",");
			}
//			bw.write("\n");
		}
//		bw.close();
	}
}
