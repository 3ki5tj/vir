#include <stdio.h>
#if defined(Macintosh) || defined(__APPLE__)
/* to compile:  gcc -framework OpenGL -framework GLUT lj3.c */
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#define ZCOM_PICK
#define ZCOM_SS
#define ZCOM_ARGOPT
#include "zcom.h"



char *fninp = NULL;
char *fnsnap = NULL, *extsnap = "png";
double *xyz; /* coordinates */
int n;
int dim, dim0;
double *dismat, *dismat0;
double *disdiff, disdiffmax = 0, disdiffmin = 0, energy = 0;

int width = 400;
int height = 400;



double modelscale = 1;

GLfloat bgcolor[] = {1.f, 1.f, 1.f, 1.f}; /* red, gree, blue, alpha */
GLfloat bglightcolor[] = {.5f, .5f, .5f, 1.f}; /* red, green, blue, alpha */
GLfloat bglightpos[] = {1.f, 2.f, 3.f, 0.f}; /* x, y, z, 0 */
int windowid = -1;
double radius = 0.5; /* radius of sphere */



static void doargs(int argc, char **argv)
{
  argopt_t *ao;

  ao = argopt_open(0);
  argopt_add(ao, NULL, "!%s", &fninp, "input position file");
  argopt_add(ao, "-w", "%d",  &width, "width of the window");
  argopt_add(ao, "-h", "%d",  &height, "height of the window");
  argopt_add(ao, "-s", "%s",  &fnsnap, "snapshot file");
  argopt_add(ao, "-x", "%s",  &extsnap, "extension of the snapshot file");
  argopt_parse(ao, argc, argv);
  if ( fninp == NULL ) argopt_help(ao);
  if ( fnsnap == NULL ) {
    char *p;
    fnsnap = ssdup(fninp);
    if ((p = strchr(fnsnap, '.')) != NULL) *p = '\0';
    sscat(fnsnap, ".");
    sscat(fnsnap, extsnap);
  }
  argopt_close(ao);
}



/* load the coordinates */
static void load(const char *fn)
{
  FILE *fp;
  char s[1024], *p;
  int ln, i = 0, j, next;

  xfopen(fp, fn, "r", return);
  if ( fgets(s, sizeof s, fp) == NULL
    || (s[0] != '#')
    || sscanf(s + 1, "%d%d%n", &dim, &n, &next) != 2 ) {
    fprintf(stderr, "%s: first line corrupted %s", fn, s);
    fclose(fp);
    return;
  }
  sscanf(s + 1 + next, "%d %lf %lf %lf",
        &dim0, &energy, &disdiffmax, &disdiffmin);
  xnew(xyz, n * dim);
  xnew(disdiff, n);
  for ( i = 0; i < n; i++ ) {
    fgets(s, sizeof s, fp);
    for ( p = s, j = 0; j < dim; j++, p += next ) {
      if ( sscanf(p, "%lf%n", &xyz[i*dim + j], &next) != 1 ) {
        fprintf(stderr, "corrupted at ln %d, i %d, j %d\n%s", ln, i, j, s);
        goto END;
      }
    }
    sscanf(p, "%lf", &disdiff[i]);
  }
  printf("file dim %d, n %d, dim0 %d, energy %g disdiffmax %g min %g\n",
      dim, n, dim0, energy, disdiffmax, disdiffmin);

  /* load the distance matrix */
  xnew(dismat, n * n);
  xnew(dismat0, n * n);
  for ( ln = 0; fgets(s, sizeof s, fp); ln++ ) {
    double dis1, dis0;
    if (s[0] != '#') break;
    sscanf(s + 1, "%d %d %lf %lf", &i, &j, &dis1, &dis0);
    dismat[i*n + j] = dis1;
    dismat0[i*n + j] = dis0;
  }
END:
  fclose(fp);
}



/* compute the scale */
static double getscale(real *x, int n, int dim)
{
  int i;
  real xm, xmax = 0;

  for ( i = 0; i < n * dim; i++ )
    if ((xm = fabs(x[i])) > xmax)
      xmax = xm;
  return xmax + .5; /* add a radus as the margin */
}



#ifdef IMLIB2

/* the following code is adapted from
 * http://abdessel.iiens.net/dev/MeshViewer/ext/MeshViewer.c */

#include <Imlib2.h>

static void gl2Imlib2(unsigned char* data, unsigned width, unsigned height)
{
  unsigned int* pixelUp, *pixelDown;
  unsigned int halfHeight = height/2;
  unsigned int i, j;
  char c;
  int pixel;
#define GL2IMLIB2(pixel) { \
  c = ((unsigned char *)(pixel))[0]; \
  ((unsigned char *)(pixel))[0] = ((unsigned char*)(pixel))[2]; \
  ((unsigned char *)(pixel))[2] = c; }

  for(i = 0; i < halfHeight; i++)
    for(j = 0; j < width; j++) {
      pixelUp = ((unsigned *) data) + (width * i + j);
      pixelDown = ((unsigned *) data) + (width * (height - 1 - i) + j);
      GL2IMLIB2(pixelUp);
      GL2IMLIB2(pixelDown);
      pixel = *pixelUp;
      *pixelUp = *pixelDown;
      *pixelDown = pixel;
    }
}

static void snapshot(const char *fn, unsigned newWidth, unsigned newHeight)
{
  GLint v[4];
  glGetIntegerv(GL_VIEWPORT, v);
  unsigned int width, height;
  unsigned char *data;
  Imlib_Image image, resizedImage;

  width = v[2] - v[0];
  height = v[3] - v[1];
  image = imlib_create_image(width, height);
  imlib_context_set_image(image);
  data = (unsigned char *) imlib_image_get_data();
  /* read pixels from the screen */
  glReadPixels(v[0], v[1], v[2], v[3], GL_RGBA, GL_UNSIGNED_BYTE, data);
  gl2Imlib2(data, width, height);
  if (newWidth <= 0) newWidth = width;
  if (newHeight <= 0) newHeight = height;
  resizedImage = imlib_create_cropped_scaled_image(
          0, 0, width, height, newWidth, newHeight);
  imlib_free_image();
  imlib_context_set_image(resizedImage);
  imlib_save_image(fn);
  imlib_free_image();
  imlib_context_set_image(image);
  fprintf(stderr, "saved the snapshot to %s\n", fn);
}
#endif /* defined(IMLIB2) */



double zoomscale = 0.85;
enum {MS_NONE, MS_ZOOM, MS_ROTATE};
int mouse_action = MS_NONE, mouse_down, mouse_x, mouse_y;



enum {
  SCALE_MINUS, SCALE_PLUS, X_MINUS, X_PLUS, Y_MINUS, Y_PLUS, Z_MINUS, Z_PLUS,
  FULLSCREEN, PRINTSCREEN, QUIT, MENULAST};

struct {
  int id;
  unsigned char key;
  const char *desc;
} menukey[MENULAST] = {
  {SCALE_MINUS, 's', "Zoom out"},
  {SCALE_PLUS,  'S', "Zoom in"},
  {X_MINUS,     'x', "Rotate around -x"},
  {X_PLUS,      'X', "Rotate around +x"},
  {Y_MINUS,     'y', "Rotate around -y"},
  {Y_PLUS,      'Y', "Rotate around +y"},
  {Z_MINUS,     'z', "Rotate around -z"},
  {Z_PLUS,      'Z', "Rotate around +z"},
  {FULLSCREEN,  'f', "Toggle full screen"},
  {PRINTSCREEN, 'p', "Save screen to file"},
  {QUIT,        'q', "Quit"}
};



static void menu(int id)
{
  GLfloat mat[4][4];

  if (id == SCALE_MINUS || id == SCALE_PLUS) {
    double s1 = zoomscale + (id == SCALE_PLUS ? 0.1 : -0.1);
    if (s1 >= 0.1) zoomscale = s1;
    glutPostRedisplay();
  } else if (id == X_MINUS || id == X_PLUS) {
    glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mat);
    glRotatef((id == X_PLUS) ? 5.f : -5.f, mat[0][0], mat[1][0], mat[2][0]);
    glutPostRedisplay();
  } else if (id == Y_MINUS || id == Y_PLUS) {
    glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mat);
    glRotatef((id == Y_PLUS) ? 5.f : -5.f, mat[0][1], mat[1][1], mat[2][1]);
    glutPostRedisplay();
  } else if (id == Z_MINUS || id == Z_PLUS) {
    glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mat);
    glRotatef((id == Z_PLUS) ? 5.f : -5.f, mat[0][2], mat[1][2], mat[2][2]);
    glutPostRedisplay();
  } else if (id == FULLSCREEN) {
    static int full = 0, x, y, w, h;
    full = !full;
    if (full) {
      x = glutGet(GLUT_WINDOW_X);
      y = glutGet(GLUT_WINDOW_Y);
      w = glutGet(GLUT_WINDOW_WIDTH);
      h = glutGet(GLUT_WINDOW_HEIGHT);
      glutFullScreen();
    } else {
      glutPositionWindow(x, y);
      glutReshapeWindow(w, h);
    }
#ifdef IMLIB2
  } else if (id == PRINTSCREEN) {
    snapshot(fnsnap, 0, 0);
#endif
  } else if (id == QUIT) {
    exit(0);
  }
}



static void initmenu(void)
{
  int i;
  char s[64];

  glutCreateMenu(menu);
  for (i = 0; i < MENULAST; i++) {
    sprintf(s, "%s, key: %c\n", menukey[i].desc, menukey[i].key);
    glutAddMenuEntry(s, menukey[i].id);
  }
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}



static void keypress(unsigned char c, int x, int y)
{
  int i;

  (void) x; (void) y;
  if (c == 27) exit(0);
  for (i = 0; i < MENULAST; i++)
    if (c == menukey[i].key)
      menu(menukey[i].id);
}



static void keyusage(void)
{
  int i;

  printf("Hotkeys:\n");
  for (i = 0; i < MENULAST; i++)
    printf(" %c: %s\n", menukey[i].key, menukey[i].desc);
}



static void display(void)
{
  float color[4] = {.2f, 1.f, 1.f, 1.f};
  double x, y, z, w, w1;
  int i;

  if ( xyz == NULL ) return;
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  for ( i = 0; i < n; i++ ) {
    float cmin = .4f;
    x = xyz[i*dim + 0];
    y = xyz[i*dim + 1];
    if (dim >= 3) {
      z = xyz[i*dim + 2];
    } else {
      z = 0;
    }

    if (disdiffmax < 1e-6) w = 0.;
    else w = (disdiff[i] - disdiffmin) / (disdiffmax - disdiffmin);

    /* set up the hue, if w falls in (0, 0.5), turn from blue to green
     * if w falls in (0.5, 1), from green to red  */
    if (w < 0.5) {
      w1 = w * 2;
      color[0] = cmin;
      color[1] = (float) (cmin + (1 - cmin) * w1);
      color[2] = (float) (1 - (1 - cmin) * w1);
    } else {
      w1 = 2*w - 1;
      color[0] = (float) (cmin + (1 - cmin) * w1);
      color[1] = (float) (1 - (1 - cmin) * w1);
      color[2] = cmin;
    }
    glLightfv(GL_LIGHT0, GL_DIFFUSE, color);
    glLightfv(GL_LIGHT0, GL_AMBIENT, color);

    glPushMatrix();
    glTranslated(x, y, z);
    glutSolidSphere(radius, 40, 40);
    glPopMatrix();
  }

  glutSwapBuffers();
}



static void reshape(int w, int h)
{
  double xs = modelscale, ys = modelscale, zs = modelscale*10;

  if (w > h) xs = modelscale * w / h;
  else ys = modelscale * h / w;
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-xs, xs, -ys, ys, -zs, zs);
  glMatrixMode(GL_MODELVIEW);
}



static void mouse(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN) {
    mouse_down++;
    if (button == 3) {
      menu(SCALE_PLUS);
    } else if (button == 4) {
      menu(SCALE_MINUS);
    }
  } else if (--mouse_down <= 0) {
    mouse_down = 0;
    mouse_action = MS_NONE;
  }
  mouse_x = x;
  mouse_y = y;
}



/* mouse motion function for GLUT */
static void motion(int x, int y)
{
  if (x == mouse_x && y == mouse_y) return;
  if (mouse_down) {
    float angx = (float)( (y - mouse_y) * 360.f / glutGet(GLUT_WINDOW_HEIGHT) );
    float angy = (float)( (x - mouse_x) * 360.f / glutGet(GLUT_WINDOW_WIDTH) );
    float mat[4][4];

    glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mat);
    glRotated(angx, mat[0][0], mat[1][0], mat[2][0]);
    glRotated(angy, mat[0][1], mat[1][1], mat[2][1]);
    glutPostRedisplay();
  }
  mouse_x = x;
  mouse_y = y;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  load(fninp);
  modelscale = getscale(xyz, n, dim) * 1.05;
  printf("model scale %g\n", modelscale);

  glutInit(&argc, argv);
  glutInitWindowSize(width, height);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  windowid = glutCreateWindow("virclus");

  /* set up the lights */
  glClearColor(bgcolor[0], bgcolor[1], bgcolor[2], bgcolor[3]);
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, bglightcolor);
  glLightfv(GL_LIGHT0, GL_POSITION, bglightpos);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);

  /* register functions */
  initmenu();
  keyusage();
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keypress);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);

  glutMainLoop();
  return 0;
}

