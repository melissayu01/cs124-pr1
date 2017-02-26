import java.util.Random;

public class EuclideanDistExp implements Experiment {
	
	Random rnd;
	
	public EuclideanDistExp() {
		this.rnd = new Random(System.currentTimeMillis());
	}
	
	@Override
	public String toString() {
		return "Euclidean Distance Experiment";
	}

	@Override
	public double weight(Vertex u, Vertex v) {
		double[] p = u.getCoords();
		double[] q = v.getCoords();			
		double weight = 0;
		for (int i = 0; i < p.length; i++) {
			weight += Math.pow(p[i] - q[i], 2);
		}
		return Math.sqrt(weight);
	}

	@Override
	public Vertex[] createGraph(int n_points, int dim) {
		Vertex[] graph = new Vertex[n_points];
		double[] coords;
		for (int i = 0; i < n_points; i++) {
			coords = new double[dim];
			for (int j = 0; j < dim; j++)
				coords[j] = rnd.nextDouble();
			graph[i] = new Vertex(i, coords);
		}
		graph[0].setMinDist(0);
		return graph;
	}
	
	public static void main(String[] args) {
		int n_points = 5;
		int dim = 3;

		EuclideanDistExp ede = new EuclideanDistExp();
		Vertex[] graph = ede.createGraph(n_points, dim);
		for (int i = 0; i < graph.length; i++) 
        	System.out.println(graph[i]);

		for (int i = 0; i < n_points; i++)
			for (int j = i+1; j < n_points; j++)
				System.out.format("Weight (%d, %d): %.3f\n", i, j, ede.weight(graph[i], graph[j]));
    }

}
