/**  
 *  Template class encapsulating core functionality for every "experiment".
 *  
 *  Classes extending this interface must implement the methods `weight(u, v)`, which 
 *  takes in 2 vertex objects and returns the edge weight between the two vertices, and
 *  `createGraph(n_points, dim)`, which creates a new random graph according to the
 *  experiment rules with `n_points` points and vertex coordinate dimension `dim` if
 *  applicable.
 *
 */
public interface Experiment {
	String toString();
	
	/**
     * Returns the weight of the edge between 2 vertices.
     *
     * @param  u one endpoint of this edge
     * @param  v other endpoint of this edge
     * @return weight of edge (u, v)
     */
	double weight(Vertex u, Vertex v);
	
	/**
     * Returns a random graph with number of points and dimension specified.
     *
     * @param  n_points number of vertices in graph
     * @param  dim dimensionality of vertex coordinates, if applicable
     * @return array of vertex objects representing a graph
     */
	Vertex[] createGraph(int n_points, int dim);
}
