import java.util.Arrays;

public class RandMST {
	
	private final Vertex[] graph;
	private int[] prioQueue;
	private int prioQueueEnd; // exclusive
	private final Experiment t;
	private final int n_trials;
	private final int n_points;	
	private final int dim;

	public RandMST(Experiment t, int n_trials, int n_points, int dim) {
		this.t = t;
		this.n_trials = n_trials;
		this.n_points = n_points;
		this.dim = dim;
		this.graph = t.createGraph(n_points, dim);
		
		this.prioQueue = new int[n_points];
		for (int i = 0; i < n_points; i++) {
			this.prioQueue[i] = i;
		}
		this.prioQueueEnd = n_points;
	}
	
	public Vertex[] getGraph() {
		return graph;
	}
	
	private Vertex vertex(int i) {
		return graph[prioQueue[i]];
	}
	
	private int parent(int i) {
		if (i < 1)
        	throw new IllegalArgumentException("vertex has no parent in heap");
		return (int) ( Math.ceil( ((double) i / 2.0) - 1) );
	}
	
	private int left(int i) {
		if (i < 0)
        	throw new IllegalArgumentException("vertex must be within valid heap boundaries");
		return 2*i + 1;
	}
	
	private int right(int i) {
		if (i < 0)
        	throw new IllegalArgumentException("vertex must be within valid heap boundaries");
		return 2*i + 2;
	}
		
	private void minHeapify(int i) {
		int l = left(i);
		int r = right(i);
		int min = i;
				
		if (l < n_points && vertex(l).compareTo(vertex(i)) < 0)
			min = l;
		if (r < n_points && vertex(r).compareTo(vertex(min)) < 0)
			min = r;
		if (min != i) {
			int tmp = prioQueue[i];
			prioQueue[i] = prioQueue[min];
			prioQueue[min] = tmp;
			minHeapify(min);
		}
	}
	
	private Vertex extractMin() {
		if (prioQueueEnd <= 0)
			throw new IllegalArgumentException("no items left in queue");
		Vertex v = vertex(0);
		v.setInMST(true);
		
		--prioQueueEnd;
		prioQueue[0] = prioQueue[prioQueueEnd];	
		minHeapify(0);
		return v;
	}
	
	private void increasePriority(int i, double newDist) {
		Vertex v = vertex(i);
		if (newDist > v.getMinDist()) 
        	throw new IllegalArgumentException("new distance must be smaller than old distance");
		v.setMinDist(newDist); 
		
		int tmp;
		while (i > 0 && vertex(parent(i)).compareTo(v) > 0) {
			tmp = prioQueue[i];
			prioQueue[i] = prioQueue[parent(i)];
			prioQueue[parent(i)] = tmp;
			i = parent(i);
			v = vertex(i);
		}
	}
	
	private double getTreeWeight() {
		double weight = 0;
		for (int i = 0; i < n_points; i++) {
			weight += graph[i].getMinDist();
		}
		return weight;
	}
	
	private void createMST() {
		Vertex u, v;
		int key;
		double newDist;
		
		System.out.println("Creating MST...");
		
		while (prioQueueEnd > 0) {
			System.out.println("\nnew iteration");
			u = extractMin();
			key = u.getKey();
			
			System.out.println(u);
			System.out.println( "current priority queue: " + Arrays.toString(
					Arrays.copyOfRange(prioQueue, 0, prioQueueEnd)) );

							
			for (int i = 0; i < prioQueueEnd; i++) {
				v = vertex(i);				
				newDist = t.weight(u, v);
				
				System.out.println("Weight (" + u.getKey() + "," + v.getKey() + "): " + newDist);
				
				if (!v.isInMST() && newDist < v.getMinDist()) {
					v.setPredecessorKey(key);
					increasePriority(i, newDist);
				}
			}
		}
		System.out.println("\nDone creating MST!");
	}
	
	public static void main(String[] args) {
		RandWeightExp exp = new RandWeightExp();
		RandMST randMST = new RandMST(exp, 1, 10, 2);
		randMST.createMST();
		double weight = randMST.getTreeWeight();
		
		Vertex[] graph = randMST.getGraph();
		for (int i = 0; i < graph.length; i++) {
			System.out.println(graph[i]);
		}
		
		System.out.println("total MST weight: " + weight);
	}

}
