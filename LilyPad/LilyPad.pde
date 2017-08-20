/*********************************************************
                  Main Window!

Click the "Run" button to Run the simulation.

Change the geometry, flow conditions, numercial parameters
visualizations and measurments from this window.

This screen has an example. Other examples are found at 
the top of each tab. Copy/paste them here to run, but you 
can only have one setup & run at a time.

*********************************************************/

Window view; // window in h-units

boolean saveimg = false;

int nx = (int)pow(2,6); // x-dir
int ny = (int)pow(2,6); // y-dir

float L = nx/6.;
float thick = 1;
float M = 10;
float stiff = 20;

float xpos = nx/4.;
float ypos = ny/2.;
PVector align = new PVector(1, 0);

float t=0;
float dt;


BDIM flow;
FlexibleSheet sheet,sheet2;
FloodPlot plot;
CollisionSolver collisions;


void setup() {
  
  size(1080, 600);
  //view = new Window(nx, ny);
  view = new Window( 1, 1, nx, ny, 0, 0, width, height);
  
  sheet = new FlexibleSheet(L, thick, M, stiff, xpos, ypos, align, view);
  
  plot = new FloodPlot(view); // standard window
  
  sheet.cpoints[0].makeFixed(); // pinning leading point
  
  dt = sheet.dtmax;
  
  flow = new BDIM(nx, ny, dt, sheet, 0.01, true);
  
  plot.range = new Scale(-1,1);
  plot.hue = new Scale(100, 40);
  plot.setLegend("pressure");
  
  //plot.range = new Scale(-1,1);
  //plot.hue = new Scale(5, 220);
  //plot.setLegend("vorticity");
  
  collisions = new CollisionSolver(sheet, view);
  
} // end of setup


void draw() {
  
  sheet.update(dt, flow);
  sheet.update2(dt, flow);
  collisions.SolveCollisions();
  
  flow.update(sheet);
  flow.update2();
  
  plot.display(flow.p);
  //plot.display(flow2.u.curl());
  
  sheet.mydisplay();
  
  t += dt;
  
  if (saveimg) saveFrame("movie/frame_######.png");
  
  //noLoop();
}