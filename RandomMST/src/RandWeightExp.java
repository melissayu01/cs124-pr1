
public class RandWeightExp implements Experiment{

	public RandWeightExp() {}

	@Override
	public String getName() {
		return "Random Weight Experiment";
	}

	@Override
	public double weight(Vertex u, Vertex v) {
		return Math.random();
	}

	@Override
	public Vertex[] createGraph(int n_points, int dim) {
		Vertex[] graph = new Vertex[n_points];
		double[] coords = {0};
		for (int i = 0; i < n_points; i++) {
			graph[i] = new Vertex(i, coords);
		}
		graph[0].setMinDist(0);
		return graph;
	}
	
	public static void main(String[] args) {
		RandWeightExp rwe = new RandWeightExp();
		Vertex[] graph = rwe.createGraph(10, 3);
		for (int i = 0; i < graph.length; i++) {
        	System.out.println(graph[i]);
		}
		for (int i = 0; i < 10; i++) {
        	System.out.format("Edge weight: %.5f\n", rwe.weight(graph[0], graph[0]));
		}
    }

}
