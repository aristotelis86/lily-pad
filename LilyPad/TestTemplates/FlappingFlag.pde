/***************************************************************
    Template to run a SINGLE flapping sheet.
    Options include:
        - defining the grid analysis nx, ny
        - Reynolds number Re
        - Mass number M
        - Length of the sheet L
        - Position of leading edge lpos
        - resolution of the sheet resol (spacing of control points)
        - stiffness of springs used
        - record the run to create a video later (saveimg)
        
***************************************************************/        
boolean saveimg = false; // save frames?

int nx = (int)pow(2,6); // x-dir
int ny = (int)pow(2,6); // y-dir

float Re = 100; // Reynolds
float L = nx/4.; // Length
float thick = 1; // thickness (display and collisions)
int resol = 1; // spacing of control points (discretization of structure)
float M = 10; // mass number (mass solid/mass fluid)
float stiff = 100; // spring stiffness
PVector lpos = new PVector( nx/3., 0.55*ny ); // position of control point 0
PVector align = new PVector(1, 0); // initial alignment of the sheet

float t=0; // time 
float dt; // time step

BDIM flow;
FlexibleSheet sheet;
CollisionHandler collide;
WriteInfo writer;
Window view; // window in h-units
FloodPlot plot;

void settings(){
    size(600, 600);
}

void setup() {
  view = new Window( 1, 1, nx, ny, 0, 0, width, height);
  
  sheet = new FlexibleSheet(L, thick, M, resol, stiff, lpos, align, view);
  sheet.cpoints[0].makeFixed(); // pinning leading point
  
  plot = new FloodPlot(view);
  
  dt = sheet.dtmax;
  
  flow = new BDIM(nx, ny, dt, sheet, 1/Re, true);
  
  plot.range = new Scale(-1,1);
  plot.hue = new Scale(100, 40);
  plot.setLegend("pressure");
  
  collide = new CollisionHandler( sheet );
  writer = new WriteInfo( sheet );
} // end of setup


void draw() {
  
  // Predict structure and fluid
  sheet.updateAlt(dt, flow);
  flow.update(sheet);
  
  // Correct structure and fluid
  sheet.updateAlt2(dt, flow);
  flow.update2();
  
  // Detect and resolve collisions
  collide.HandleCollisions();
  
  // Display
  plot.display(flow.p);
  sheet.display();
  
  // Create output files (txt and png)
  writer.dampAllInfo( t, false );
  if (saveimg) saveFrame("movie/frame_######.png");
  
  t += dt;
}

// Ensure that the txt files are closed properly
void keyPressed() {
  writer.terminateFiles();
}


