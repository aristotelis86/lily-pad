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
FloodPlot plot;

boolean saveimg = false;

int nx = (int)pow(2,6); // x-dir
int ny = (int)pow(2,6); // y-dir

float Re = 100;
float L = nx/8.;
float thick = 1;
int resol = 1;
float M = 10;
float stiff = 100;
PVector lpos = new PVector( nx/3., 0.55*ny );
PVector align = new PVector(1, 0);

float t=0;
float dt;

BDIM flow;
FlexibleSheet sheet;
CollisionHandler collide;
WriteInfo writer;

void settings(){
    size(600, 600);
}

void setup() {
  view = new Window( 1, 1, nx, ny, 0, 0, width, height);
  
  sheet = new FlexibleSheet(L, thick, M, resol, stiff, lpos, align, view);
  sheet.cpoints[0].makeFixed(); // pinning leading point
  
  plot = new FloodPlot(view); // standard window
  
  dt = sheet.dtmax;
  
  flow = new BDIM(nx, ny, dt, sheet, sheet.Length/Re, true);
  
  plot.range = new Scale(-1,1);
  plot.hue = new Scale(100, 40);
  plot.setLegend("pressure");
  
  collide = new CollisionHandler( sheet );
  writer = new WriteInfo( sheet );
} // end of setup


void draw() {
  
  sheet.updateAlt(dt, flow);
  flow.update(sheet);
  
  sheet.updateAlt2(dt, flow);
  flow.update2();
  
  collide.HandleCollisions();
  
  plot.display(flow.p);
  sheet.display();
  
  writer.dampAllInfo( t, false );
  if (saveimg) saveFrame("movie/frame_######.png");
  t += dt;
}

void keyPressed() {
  writer.terminateFiles();
}