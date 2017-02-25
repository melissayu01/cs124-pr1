
public interface Experiment {
	String toString();
	double weight(Vertex u, Vertex v);
	Vertex[] createGraph(int n_points, int dim);
}
