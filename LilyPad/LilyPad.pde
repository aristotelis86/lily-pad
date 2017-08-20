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

float L = nx/8.;
float thick = 1;
float M = 10;
float stiff = 20;

float xpos = nx/6.;
float xpos2 = 0.6*nx;
float ypos = 20;
float ypos2 = 25;
PVector align = new PVector(1, 0);

float t=0;
float dt;


BDIM flow;
BodyUnion bodies;
FlexibleSheet [] sheet;
Body circle;
CollisionSolver collisions;


void setup() {
  
  //size(1080, 600);
  size(600, 600);
  //view = new Window(nx, ny);
  view = new Window( 1, 1, nx, ny, 0, 0, width, height);
  
  sheet = new FlexibleSheet[4];
  
  
  sheet[0] = new FlexibleSheet(L, thick, M, stiff, xpos, ypos, align, view);
  sheet[1] = new FlexibleSheet(L, thick, M, stiff, xpos, ypos2, align, view);
  sheet[2] = new FlexibleSheet(L, thick, M, stiff, xpos2, ypos, align, view);
  sheet[3] = new FlexibleSheet(L, thick, M, stiff, xpos2, ypos2, align, view);
  
  circle = new CircleBody(nx/3. + nx/12., ny/3., nx/6., view);
  
  plot = new FloodPlot(view); // standard window
  
  sheet[0].cpoints[0].makeFixed(); // pinning leading point
  sheet[1].cpoints[0].makeFixed(); // pinning leading point
  sheet[2].cpoints[0].makeFixed(); // pinning leading point
  sheet[3].cpoints[0].makeFixed(); // pinning leading point
  
  dt = sheet[0].dtmax;
  
  bodies = new BodyUnion(sheet[0], sheet[1]);
  bodies.add(sheet[2]);
  bodies.add(sheet[3]);
  bodies.add(circle);
  
  
  flow = new BDIM(nx, ny, dt, bodies, 0.01, true);
  
  plot.range = new Scale(-1,1);
  plot.hue = new Scale(100, 40);
  plot.setLegend("pressure");
  
  //plot.range = new Scale(-1,1);
  //plot.hue = new Scale(5, 220);
  //plot.setLegend("vorticity");
  
  collisions = new CollisionSolver(sheet, view);
  
} // end of setup


void draw() {
  
  for (FlexibleSheet fl : sheet) {
    fl.update(dt, flow);
    fl.update2(dt, flow);
  }
  
  collisions.SolveCollisions();
  
  flow.update(bodies);
  flow.update2();
  
  plot.display(flow.p);
  
  
  for (FlexibleSheet fl : sheet) fl.display();
  circle.display();
  
  t += dt;
  
  if (saveimg) saveFrame("movie/frame_######.png");
  
  //noLoop();
}