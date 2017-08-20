/*********************************************************
                  Main Window!

Click the "Run" button to Run the simulation.

Change the geometry, flow conditions, numercial parameters
visualizations and measurments from this window.

This screen has an example. Other examples are found at 
the top of each tab. Copy/paste them here to run, but you 
can only have one setup & run at a time.

*********************************************************/

Window view, view2; // window in h-units
FloodPlot plot, plot2;

boolean saveimg = true;

int nx = (int)pow(2,6); // x-dir
int ny = (int)pow(2,6); // y-dir

float L = nx/6.;
float thick = 1;
float M = 10;
float stiff = 20;

float xpos = nx/4.;
float ypos = 20;
float ypos2 = 25;
PVector align = new PVector(1, 0);

float t=0;
float dt;


BDIM flow, flow2;
BodyUnion bodies, bodies2;
FlexibleSheet [] sheet, sheet2;
CollisionSolver collisions, collisions2;


void setup() {
  
  size(1080, 600);
  //view = new Window(nx, ny);
  view = new Window( 1, 1, nx, ny, 0, 0, width, height/2);
  view2 = new Window( 1, 1, nx, ny, 0, height/2, width, height/2);
  
  sheet = new FlexibleSheet[2];
  sheet2 = new FlexibleSheet[2];
  
  sheet[0] = new FlexibleSheet(L, thick, M, stiff, xpos, ypos, align, view);
  sheet[1] = new FlexibleSheet(L, thick, M, stiff, xpos, ypos2, align, view);
  sheet2[0] = new FlexibleSheet(L, thick, M, stiff, xpos, ypos, align, view2);
  sheet2[1] = new FlexibleSheet(L, thick, M, stiff, xpos, ypos2, align, view2);
  
  plot = new FloodPlot(view); // standard window
  plot2 = new FloodPlot(view2); // standard window
  
  //sheet[0].cpoints[0].makeFixed(); // pinning leading point
  //sheet[1].cpoints[0].makeFixed(); // pinning leading point
  //sheet2[0].cpoints[0].makeFixed(); // pinning leading point
  //sheet2[1].cpoints[0].makeFixed(); // pinning leading point
  
  dt = sheet[0].dtmax;
  
  bodies = new BodyUnion(sheet[0], sheet[1]);
  bodies2 = new BodyUnion(sheet2[0], sheet2[1]);
  
  flow = new BDIM(nx, ny, dt, bodies, 0.01, true);
  flow2 = new BDIM(nx, ny, dt, bodies2, 0.01, true);
  
  plot.range = new Scale(-1,1);
  plot.hue = new Scale(100, 40);
  plot.setLegend("pressure");
  
  plot2.range = new Scale(-1,1);
  plot2.hue = new Scale(5, 220);
  plot2.setLegend("vorticity");
  
  collisions = new CollisionSolver(sheet, view);
  collisions2 = new CollisionSolver(sheet2, view2);
  
} // end of setup


void draw() {
  
  for (FlexibleSheet fl : sheet) {
    fl.update(dt, flow);
    fl.update2(dt, flow);
  }
  for (FlexibleSheet fl : sheet2) {
    fl.update(dt, flow2);
    fl.update2(dt, flow2);
  }
  
  collisions.SolveCollisions();
  collisions2.SolveCollisions();
  
  flow.update(bodies);
  flow.update2();
  
  flow2.update(bodies2);
  flow2.update2();
  
  plot.display(flow.p);
  plot2.display(flow.u.curl());
  
  for (FlexibleSheet fl : sheet) fl.display();
  for (FlexibleSheet fl : sheet2) fl.display();
  
  t += dt;
  
  if (saveimg) saveFrame("movie/frame_######.png");
  
  //noLoop();
}