/***************************************************************
    Template to run TWO flapping sheets.
    Options include:
        - defining the grid analysis nx, ny
        - Reynolds number Re
        - Mass number M
        - Length of the sheet L
        - Position of the centre between the 
          leading edges of the two sheets
        - Separation between the sheets separ  
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
int resol = 1; // discretization of structures
float M = 5; // Mass number (mass solid/mass fluid)
float stiff = 100; // stiffness of spring
PVector lpos = new PVector( nx/3., 0.5*ny ); // center position
PVector align = new PVector(1, 0); // initial alignment of the sheets
float separ = L/2.; // separation of the sheets

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
  sheet = new FlexibleSheet[2];
  
  sheet[0] = new FlexibleSheet(L, thick, M, resol, stiff, new PVector(lpos.x, lpos.y-separ/2.), align, view);
  sheet[0].cpoints[0].makeFixed(); // pinning leading point
  sheet[1] = new FlexibleSheet(L, thick, M, resol, stiff, new PVector(lpos.x, lpos.y+separ/2.), align, view);
  sheet[1].cpoints[0].makeFixed(); // pinning leading point
  
  plot = new FloodPlot(view); // standard window
  
  dt = sheet[0].dtmax;
  bodies = new BodyUnion(sheet[0], sheet[1]);
  
  flow = new BDIM(nx, ny, dt, bodies, 1/Re, true);
  
  plot.range = new Scale(-1,1);
  plot.hue = new Scale(100, 40);
  plot.setLegend("pressure");
  
  collide = new CollisionHandler( sheet );
  writer = new WriteInfo( sheet );
} // end of setup


void draw() {
  
  // Predict structure and fluid
  for (FlexibleSheet fs : sheet) fs.updateAlt(dt, flow);
  flow.update(bodies);
  
  // Correct structure and fluid
  for (FlexibleSheet fs : sheet) fs.updateAlt2(dt, flow);
  flow.update2();
  
  // Detect and resolve collisions
  collide.HandleCollisions();
  
  // Display
  plot.display(flow.p);
  for (FlexibleSheet fs : sheet) fs.display();
  
  // Produce output (txt and png)
  writer.dampAllInfo( t, false );
  if (saveimg) saveFrame("movie/frame_######.png");
  
  t += dt;
}

// Ensure proper closing of txt files
void keyPressed() {
  writer.terminateFiles();
}


