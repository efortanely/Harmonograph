/*
////////////////////////////////
Right-button/Space will generate a new lissajous curve
Left-button will replay the current lissajous curve
////////////////////////////////
*/

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.event.KeyEvent;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.util.Arrays;
import java.util.Random;

public class FunctionDisplayer extends Engine{
	Harmonograph harmonograph;
	Path2D.Double points;
	int maxTime = 200;
	double time;
	double timeDif = .01;
	
	@Override
	public void first() {
		setFps(500);
		initializeHarmonograph();
	}

	Point2D.Double prevPt;
	@Override
	public void tick() {
		if(newHarmonograph){
			initializeHarmonograph();
			newHarmonograph = false;
		}else if(resetHarmonograph){
			points.reset();
			time = 0;
			resetHarmonograph = false;
		}else{
			if(time < maxTime){
				Point2D.Double curPt = new Point2D.Double(harmonograph.getX(time), harmonograph.getY(time));
				if(time==0){ 
					points.moveTo(curPt.x, curPt.y);
				}else{
					points.quadTo(prevPt.x, prevPt.y, curPt.x, curPt.y);
				}
				prevPt = curPt;
			}
			time += timeDif;
		}
	}

	@Override
	public void render(Graphics2D g) {
		setSmooth(true);
		g.setColor(new Color(0,0,0,220));
		g.fillRect(0, 0, getWidth(), getHeight());
		
		g.translate(getWidth()/2, getHeight()/2);
		Color transBlue = new Color(30,50,255,150);
		g.setColor(transBlue);
		g.draw(points);
	}

	boolean newHarmonograph = false;
	boolean resetHarmonograph = false;
	@Override public void keyPressed(KeyEvent e){
		switch(e.getKeyCode()){
			case KeyEvent.VK_SPACE:
				newHarmonograph = true;
				break;
			case KeyEvent.VK_LEFT:
				resetHarmonograph = true;
				break;
			case KeyEvent.VK_RIGHT:
				newHarmonograph = true;
			break;
				
		}
	}
	
	public void initializeHarmonograph(){
		points = new Path2D.Double();
		harmonograph = new Harmonograph();
		time = 0;
	}
	
	public FunctionDisplayer(int width, int height, String name) {
		super(width,height,name);
	}
	
	public static void main(String[] s){
		new FunctionDisplayer(500,500,"Harmonograph").start();
	}
	
	private class Harmonograph{
		int[] amplitude = new int[4];
		int[] frequency = new int[4];
		double[] phase = new double[4];
		double[] damping = new double[4];
		boolean userDefinedInputs = !true;
		
		public Harmonograph(){
			if(userDefinedInputs){
				/*
				//My anxiety:
				amplitude = new int[]{60,60,120,60};
				frequency = new int[]{6,9,9,5};
				phase = new double[]{1.6060233193580584,4.8307845414828705,6.039077722161357,0.9579844763894997};
				damping = new double[]{0.046765235202287724,0.008663245926084057,0.017405749831702733,0.030522271792354304};
				*/
				/*
				//Vitruvian Man:
				amplitude = new int[]{60,120,60,120};
				frequency = new int[]{9,3,6,4};
				phase = new double[]{0.2779006829265938,5.621794251322227,4.833030205245015,1.9460906245557175};
				damping = new double[]{0.021495953226130664,0.043975750204710604,0.047784599849549365,0.020607864196932085};
				*/
				/*
				//Wonder Woman symbol and/or bug:
				amplitude = new int[]{120,60,120,60};
				frequency = new int[]{1,9,6,8};
				phase = new double[]{4.013186607907121,1.3309670325473504,1.599151007470414,2.9279279101473765};
				damping = new double[]{0.009711502930995968,0.0158714159140697,0.016977220472395454,0.020679937759828425};
				*/
			}else{
				generateRandomHarmonograph();	
			}
		}
		
		public void generateRandomHarmonograph(){
			Random rand = new Random();
			for(int i = 0; i < 4; i++){
				amplitude[i] = 25*(rand.nextInt(4)+2);
				frequency[i] = 3*(rand.nextInt(4)+1);
				phase[i] = rand.nextDouble()*2*Math.PI;
				damping[i] = rand.nextDouble() / 20.0;
			}
			
			System.out.println("New harmonograph inputs!");
			String amplitudeOutput = Arrays.toString(amplitude);
			System.out.println("amplitude = new int[]{" + amplitudeOutput.substring(1, amplitudeOutput.length()-1) + "};");
			String frequencyOutput = Arrays.toString(frequency);
			System.out.println("frequency = new int[]{" + frequencyOutput.substring(1, frequencyOutput.length()-1) + "};");
			String phaseOutput = Arrays.toString(phase);
			System.out.println("phase = new double[]{" + phaseOutput.substring(1, phaseOutput.length()-1) + "};");
			String dampingOutput = Arrays.toString(damping);
			System.out.println("damping = new double[]{" + dampingOutput.substring(1, dampingOutput.length()-1) + "};");
			System.out.println();
		}
		
		public double getX(double time){
			return component(time,0) + component(time,1);
		}
		
		public double getY(double time){
			return component(time,2) + component(time,3);
		}
		
		private double component(double time, int index){
			return amplitude[index] * Math.sin(time * frequency[index] + 
				   phase[index]) * Math.exp(-damping[index] * time);
		}
	}
}