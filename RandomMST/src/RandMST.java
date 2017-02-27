import java.util.Arrays;

public class RandMST {
	
	private final Vertex[] graph;
	private int[] prioQueue;
	private int prioQueueEnd; // exclusive
	private final Experiment t;
	private final int n_points;	
	private final int dim;
	private final int flag;

	public RandMST(int n_points, int dim, int flag) {
		this.n_points = n_points;
		this.dim = dim;
		this.flag = flag;
		
		if (dim == 0) this.t = new RandDistExp();
		else this.t = new EuclideanDistExp();
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
				
		if (l < prioQueueEnd && vertex(l).compareTo(vertex(i)) < 0)
			min = l;
		if (r < prioQueueEnd && vertex(r).compareTo(vertex(min)) < 0)
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
		
		prioQueueEnd -= 1;
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
		for (int i = 0; i < n_points; i++)
			weight += graph[i].getMinDist();
		return weight;
	}
	
	private void createMST() {
		Vertex u, v;
		int key;
		double newDist;
		long startTime = System.currentTimeMillis();
		
		if (flag > 0)
			System.out.println("Creating MST...");
		
		while (prioQueueEnd > 0) {
			if (flag > 1 && (n_points - prioQueueEnd) % 10000 == 0)
				System.out.format("\n%d / %d vertices added to MST\n", n_points - prioQueueEnd, n_points);
			u = extractMin();
			key = u.getKey();
			
			if (flag > 1) {
				System.out.format("\n%s\n", u);
				System.out.println( "Current priority queue: " + Arrays.toString(
						Arrays.copyOfRange(prioQueue, 0, prioQueueEnd)) );
			}
							
			for (int i = 0; i < prioQueueEnd; i++) {
				v = vertex(i);				
				newDist = t.weight(u, v);
				
				if (flag > 1) 
					System.out.format("Weight (%d, %d): %.4f\n", u.getKey(), v.getKey(), newDist);
				
				if (newDist < v.getMinDist()) {
					v.setPredecessorKey(key);
					increasePriority(i, newDist);
				}
			}
		}
		long time = System.currentTimeMillis() - startTime;
		if (flag > 0) 
			System.out.format("...Finished creating MST in %d seconds\n", time / 1000);
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < graph.length; i++) {
			sb.append(graph[i]);
			if (i != graph.length-1)
				sb.append("\n");
		}
		return sb.toString();
    }
	
	public static double runExperiment(int flag, int n_points, int n_trials, int dim) {
		long startTime = System.currentTimeMillis();
		double totalWeight = 0;
		double weight;
		RandMST randMST;
		for (int i = 0; i < n_trials; i++) {
			randMST = new RandMST(n_points, dim, flag);
			randMST.createMST();
			weight = randMST.getTreeWeight();
			totalWeight += weight;
			
			if (flag > 1) {
				System.out.println("\nFinal MST:");
				System.out.println(randMST);
			}
			if (flag > 0)
				System.out.format("Current MST weight: %.4f\n\n", weight);
		}
		System.out.format("=========== %d MSTs created in %d minutes ===========\n", n_trials, (System.currentTimeMillis() - startTime) / (60 * 1000));
		System.out.format("%.4f %d %d %d\n\n", 
				totalWeight / n_trials, n_points, n_trials, dim);
		return totalWeight / n_trials;
	}
	
	public static void main(String[] args) {
		int flag = 0;
		int n_points = 0;
		int n_trials = 0;
		int dim = 0;
		
		if (args.length == 4) {
		    try {
		    	flag = Integer.parseInt(args[0]);
		    	n_points = Integer.parseInt(args[1]);
		    	n_trials = Integer.parseInt(args[2]);
		    	dim = Integer.parseInt(args[3]);
		    } catch (NumberFormatException ex) {
		        System.err.println("Arguments must be integers");
		        System.exit(1);
		    }
		} else {
			System.err.println("Usage: ./randmst flag n_points n_trials dim");
			System.exit(1);
		}
		
		runExperiment(flag, n_points, n_trials, dim);
	}
}
