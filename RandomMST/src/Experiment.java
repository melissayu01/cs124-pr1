
public interface Experiment {
	String getName();
	double weight(Vertex u, Vertex v);
	Vertex[] createGraph(int n_points, int dim);
}
