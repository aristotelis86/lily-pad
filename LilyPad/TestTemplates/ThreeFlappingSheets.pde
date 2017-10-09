/***************************************************************
    Template to run THREE flapping sheets.
    Options include:
        - defining the grid analysis nx, ny
        - Reynolds number Re
        - Mass number M
        - Length of the sheet L
        - Position of the leading edge of the middle sheet
        - Separation between the sheets separ  
        - resolution of the sheet resol (spacing of control points)
        - stiffness of springs used
        - record the run to create a video later (saveimg)
        
***************************************************************/

boolean saveimg = false;

int nx = (int)pow(2,6); // x-dir
int ny = (int)pow(2,6); // y-dir

float Re = 100; // Reynolds
float L = nx/4.; // Length
float thick = 1; // thickness (display and collisions)
int resol = 1; // discretization of structures
float M = 5; // mass number (mass solid/mass fluid)
float stiff = 100; // stiffness of spring
PVector lpos = new PVector( nx/3., 0.5*ny ); // centre position
PVector align = new PVector(1, 0); // initial alignment of sheets
float separ = L/2.; // separation

float t=0; // time
float dt; // time step

BDIM flow;
FlexibleSheet [] sheet;
BodyUnion bodies;
CollisionHandler collide;
WriteInfo writer;
Window view; // window in h-units
FloodPlot plot;

void settings(){
    size(600, 600);
}

void setup() {
  view = new Window( 1, 1, nx, ny, 0, 0, width, height);
  sheet = new FlexibleSheet[3];
  
  sheet[0] = new FlexibleSheet(L, thick, M, resol, stiff, new PVector(lpos.x, lpos.y-separ), align, view);
  sheet[0].cpoints[0].makeFixed(); // pinning leading point
  sheet[1] = new FlexibleSheet(L, thick, M, resol, stiff, new PVector(lpos.x, lpos.y), align, view);
  sheet[1].cpoints[0].makeFixed(); // pinning leading point
  sheet[2] = new FlexibleSheet(L, thick, M, resol, stiff, new PVector(lpos.x, lpos.y+separ), align, view);
  sheet[2].cpoints[0].makeFixed(); // pinning leading point
  
  plot = new FloodPlot(view);
  
  dt = sheet[0].dtmax;
  bodies = new BodyUnion(sheet[0], sheet[1]);
  bodies.add(sheet[2]);
  
  flow = new BDIM(nx, ny, dt, bodies, 1/Re, true);
  
  plot.range = new Scale(-1,1);
  plot.hue = new Scale(100, 40);
  plot.setLegend("pressure");
  
  collide = new CollisionHandler( sheet );
  writer = new WriteInfo( sheet );
} // end of setup


void draw() {
  
  // Predict structures and fluid
  for (FlexibleSheet fs : sheet) fs.updateAlt(dt, flow);
  flow.update(bodies);
  
  // Correct structures and fluid
  for (FlexibleSheet fs : sheet) fs.updateAlt2(dt, flow);
  flow.update2();
  
  // Detect and resolve collisions
  collide.HandleCollisions();
  
  // Display
  plot.display(flow.p);
  for (FlexibleSheet fs : sheet) fs.display();
  
  // Create output files (txt and png)
  writer.dampAllInfo( t, false );
  if (saveimg) saveFrame("movie/frame_######.png");
  
  t += dt;
}

// Ensure proper closing of txt files
void keyPressed() {
  writer.terminateFiles();
}


