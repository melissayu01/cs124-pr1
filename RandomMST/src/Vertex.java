import java.util.Arrays;

public class Vertex implements Comparable<Vertex>{
	
    private double minDist = Double.POSITIVE_INFINITY;
    private int predecessorKey = -1;
    private boolean isInMST = false;
    private final double[] coords;
    private final int key;
    
    public Vertex(int key, double[] coords) {
        if (key < 0) 
        	throw new IllegalArgumentException("vertex key must be a nonnegative integer");
        this.key = key;
        this.coords = coords;
    }

	public double getMinDist() {
		return minDist;
	}

	public void setMinDist(double minDist) {
		this.minDist = minDist;
	}

	public int getPredecessorKey() {
		return predecessorKey;
	}

	public void setPredecessorKey(int predecessorKey) {
		this.predecessorKey = predecessorKey;
	}

	public boolean isInMST() {
		return isInMST;
	}

	public void setInMST(boolean isInMST) {
		this.isInMST = isInMST;
	}
	
	public int getKey() {
		return key;
	}
	
	@Override
    public int compareTo(Vertex that) {
        return Double.compare(this.minDist, that.minDist);
    }

    public String toString() {
        return String.format("V%d %s: (%.5f) (%d=>)", key, Arrays.toString(coords), 
        		minDist, predecessorKey);
    }

    public static void main(String[] args) {
    	double[] coords = {0.2, 0.66};
        Vertex v = new Vertex(1, coords);
        v.setMinDist(1);
        Vertex u = new Vertex(2, coords);
        u.setMinDist(0.5);
        System.out.println("u = " + u);
        System.out.println("v = " + v);
        System.out.println("u < v: " + (u.compareTo(v) == -1));
    }

}
