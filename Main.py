from pygame.locals import DOUBLEBUF #pygl
from pygame import quit, QUIT, KEYUP, K_SPACE, init, OPENGL
from pygame.display import set_mode, flip
from pygame.font import SysFont
from pygame.event import get
from pygame.image import tostring
from pygame.time import wait

from OpenGL.GL import glClearColor, glEnable, GL_DEPTH_TEST,\
     GL_BLEND, glBlendFunc, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA,\
     glClearDepth, glDepthFunc, GL_LESS, glShadeModel, GL_SMOOTH,\
     glTranslatef, glClear, GL_COLOR_BUFFER_BIT,\
     GL_DEPTH_BUFFER_BIT, glRotate, glPushMatrix, glTranslate,\
     glRotatef, glBegin, GL_TRIANGLE_FAN, glColor, glVertex,\
     glEnd, GL_TRIANGLE_STRIP, glPopMatrix, glLineWidth,\
     GL_LINES, glLoadIdentity, glColor3f, glVertex3f,\
     glRasterPos2d, glDrawPixels, GL_RGBA, GL_UNSIGNED_BYTE,\
     GL_LINE_STRIP, glColor4f, GL_POLYGON
from OpenGL.GLU import gluPerspective

from utils.constants import WINDOW_SIZE, NUM_FRAMES
from numpy import pi, float64, zeros, array

from utils.constants import GRAPH_INTERVAL, DEFAULT_DELAY,\
     GREEN3F, BLUE3F, WHITE4F, GRAPH_POINTS, WINDOW_SIZE,\
     YELLOW3F, GREEN3F, BLACK3F

from math import pi, cos, sin

x, y, z = (1, 0, 0), (0, 1, 0), (0, 0, 1)
a_X, a_Y, a_Z = array(x), array(y), array(z)

class Axes:
    def __init__(self, pos):
        self.pos, (self.x, self.y, self.z) = array(pos),\
                          (3 * a_X, 3 * a_Y, 3 * a_Z)

    """
    Axes._draw_line(vector:np.array, color:tuple)
        draws a vector from its position
    """
    def _draw_line(self, vector, color):
        glColor(color)
        glBegin(GL_LINES)

        (*(glVertex(self.pos + v) for v in (zeros(3), vector)),)

        glEnd()

    """
    Axes.render([pos:tuple])
        renders the x, y and z axes
    """
    def render(self, pos=zeros(3)):     
        glPushMatrix()

        self.pos = array(pos)

        glLineWidth(2)

        #ndarray no hasheable
        (*(self._draw_line(eje, color) for eje, color in
            ((self.x, YELLOW3F), (self.y, GREEN3F),
             (self.z, BLUE3F))),)

        glLineWidth(1)
        glPopMatrix()

def draw_cylinder1(v, v1, circle_pts, z1):   
    glBegin(GL_TRIANGLE_FAN)
    glColor(*v)

    ##################
    glVertex(*(v1[2]*a_Z))

    (*(glVertex(*v, z1) for v in circle_pts),)       

    glEnd()

X, Y, Z = (1,0,0), (0,1,0), (0,0,1)

def draw_cylinder(radius, height, num_slices):
    r, n, z, circle_pts = radius, float(num_slices), height / 2,\
                          []
    
    for i in range(int(n) + 1):
        angle = 2 * pi * i/n
        x, y = r * cos(angle), r * sin(angle)
        pt = (x, y)
        circle_pts.append(pt)

    (*(draw_cylinder1(eje, (x, y, z), circle_pts, n) for eje, n
       in ((X,z), (Z,-z))),)
    
    glBegin(GL_TRIANGLE_STRIP)#draw the tube
    glColor(*Y)

    '''v = x,y

    (*((*(glVertex(*v, z) for v, z in ((v, z), (v,-z))),) for v in circle_pts),)'''
    
    (*((*(glVertex(x, y, z) for z in (z, -z)),) for x, y in circle_pts),)
        
    glEnd()

class Cylinder:
    def __init__(self, pos, rotation, height, radius, mass):
        self.rotation, self.pos, self.height, self.radius,\
            self.mass, self.slices = array(rotation), array(pos),\
            height, radius, mass, 10

    """
    Cylinder.render()
        renders the cylinder
    """
    def render(self):
        glPushMatrix()
        glTranslate(*self.pos)
     
        (*(glRotatef(self.rotation[i], *(x, y, z)[i])
           for i in range(3)),)
        
        draw_cylinder(self.radius, self.height, self.slices)

        glPopMatrix()

class Tbar:
    def __init__(self, size=2.5, pos=zeros(3),
                 angvel=array((40, 1, 1))/10):
        self.size = size

        #no colapsable
        self.initial_angvel = array(angvel)
        self.angvel, self.angvels, self.angacc, self.pos =\
            self.initial_angvel.copy(),\
            [self.initial_angvel.copy()], zeros(3), array(pos)

        # positions and rotations relative to the tbar

        self.handle, self.axis = Cylinder(rotation=90*array(x),
            pos=zeros(3), radius=self.size/3,
            height=2 * self.size, mass=4/3 * self.size),\
            Cylinder(rotation=90*a_Y, pos=self.size / 2 * a_X,
            radius=self.size/4, height=self.size,
            mass=self.size/3)

        self.cm = array((-(self.axis.height + 2 *\
                           self.handle.radius)/5, 0, 0))
        self.axes, self.graph_count, self.angvels_count =\
                   Axes(pos=self.cm), 0, 1

        self._compute_moment_inertia()

    """
    Tbar.step()
        computes the new velocities
    """   
    def step(self):
        self._compute_angacc()
        self._compute_angvel()
        self._compute_rotation()

    """
    Tbar.render()
        renders all the components (handle, axis, and axes)
    """
    def render(self):
        self.handle.render()
        self.axis.render()
        self.axes.render()

    """
    Tbar._compute_angacc()
        computes the new angular acceleration based on Euler's equations
    """
    def _compute_angacc(self):
        self.angacc = (((self.moment_inertia[1] - self.moment_inertia[2])\
            / self.moment_inertia[0] * self.angvel[1] *
            self.angvel[2] ),\
            ((self.moment_inertia[2] - self.moment_inertia[0])\
            / self.moment_inertia[1] * self.angvel[2] *
            self.angvel[0] ),\
            ((self.moment_inertia[0] - self.moment_inertia[1])\
            / self.moment_inertia[2] * self.angvel[0] *
            self.angvel[1] ))
         
    """
    Tbar._compute_angvel()
        computes the new angular velocity based on the angular acceleration
        and appends it to the plot
    """
    def _compute_angvel(self):
        for axis in range(3): self.angvel[axis] +=\
            self.angacc[axis]/DEFAULT_DELAY

        if self.graph_count == GRAPH_INTERVAL:
            self.graph_count = 0
            self.angvels_count += 1

            self.angvels.append(self.angvel.copy())

            if self.angvels_count >= GRAPH_POINTS: del\
               self.angvels[0] #OJO!

        self.graph_count += 1

    """
    Tbar._compute_rotation()
        rotates the tbar based on the angular velocity
    """
    def _compute_rotation(self):
        angx, angy, angz = self.angvel / DEFAULT_DELAY * 180 / pi

        (*(glRotate(ang, *eje) for ang, eje in ((angx, x),\
            (angy, y), (angz, z))),)

    """
    Tbar._compute_moment_inertia()
        computes the moment of inertia of the tbar
    """
    def _compute_moment_inertia(self):
        '''
        h -> v
        0 -> 2
        1 -> 0
        2 -> 1

        distance between center of mass of the handle and the total object
        '''
        self.cm_distance = self.cm[0]
        self.moment_inertia = array((self.axis.mass/2 * self.axis.radius ** 2\
               + self.handle.mass/4 * self.handle.radius ** 2\
               + self.handle.mass/12 * self.handle.height ** 2,\
               self.handle.mass/2 * self.handle.radius ** 2\
               + self.handle.mass * self.cm_distance ** 2 +\
               self.axis.mass / 4 * self.axis.radius ** 2\
               + self.axis.mass / 12 * self.axis.height ** 2\
               + self.axis.mass * (self.axis.height / 2 +\
               self.handle.radius - self.cm_distance),\
               self.handle.mass/4 * self.handle.radius ** 2\
               + self.handle.mass/12 * self.handle.height ** 2\
               + self.handle.mass * self.cm_distance ** 2\
               + self.axis.mass/4 * self.axis.radius ** 2\
               + self.axis.mass/12 * self.axis.height ** 2\
               + self.axis.mass * (self.axis.height / 2 +\
               self.handle.radius - self.cm_distance)))

class World:
    def __init__(self, window_size=WINDOW_SIZE, clear_color=WHITE4F):
        (self.width, self.height), self.clear_color,\
                     self.objects = window_size, clear_color, [] 
        self.perspective = (60, self.width / self.height, 0.1,\
                            100)

    """
    World.init_scene()
        sets some parameters on opengl
    """
    def init_scene(self):
        glClearColor(*self.clear_color)

        (*(glEnable(param) for param in (GL_DEPTH_TEST, GL_BLEND)),)
        
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

        glClearDepth(1)
        glDepthFunc(GL_LESS)
        glShadeModel(GL_SMOOTH)

        gluPerspective(*self.perspective)
        glTranslatef(*(-10*a_Z))

    """
    World.add_object(obj:any)
        adds an object to the rendering list
    """
    def add_object(self, obj):    
        self.objects.append(obj)

    """
    World.step()
        steps all the objects on the rendering list
    """
    def step(self):
        (*(obj.step() for obj in self.objects),)
            

    """
    World.render()
        renders all the objects on the rendering list
    """
    def render(self):
        (*(obj.render() for obj in self.objects),)
            

    """
    World.clear()
        clears the color and depth buffer
    """
    def clear(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

class Screen:
    def __init__(self, display, tbar, window_size=WINDOW_SIZE):
        # display
        self.display, self.window_size = display, window_size

        # menu
        self.menu_h, self.menu_padding = 1, 0.2

        # font
        self.font = SysFont("FreeMono", 16)

        # graph
        self.colors, self.graph_ratio =\
            (YELLOW3F, GREEN3F, BLUE3F),\
            ((self.menu_h/2 - self.menu_padding)/
            (tbar.initial_angvel[0] * tbar.moment_inertia[0] + 10))
        self.time_interval, self.tbar = 2 / GRAPH_POINTS *\
                    (1 - self.menu_padding), tbar

    """
    Screen.render()
        renders both the graph and the black overlay
    """
    def render(self):
        glPushMatrix()
        glLoadIdentity()

        self._render_graph()
        self._render_menu()

        glPopMatrix()

    """
    Screen.show_paused_message()
        displays a paused message on the screen
    """
    def show_paused_message(self):
        glPushMatrix()
        glLoadIdentity()
        self._render_text(array((8, 9))/10, "Paused")
        glPopMatrix()

        flip()

    """
    Screen._render_menu()
        renders the black overlay using the properties menu_h and menu_padding
    """
    def _render_menu(self):
        glColor4f(*BLACK3F, 0.25)

        glBegin(GL_POLYGON)

        (*(glVertex3f(*v) for v in (-array((1, 1, 0)),
            (1, -1, 0), array((1, self.menu_h, 0))-a_Y,
            array((-1, self.menu_h, 0))-a_Y)),)
        
        glEnd()

    """
    Screen._render_graph()
        plots the tbar's angular momentum
    """
    def _render_graph(self):
        glColor3f(*BLACK3F)

        glBegin(GL_LINES)
        # vertical line
        M, N = self.menu_padding-1, -0.5

        (*(glVertex3f(m, n, 0) for m, n in
           ((M, self.menu_padding-1.1),(M, 0.1-self.menu_padding),

            (self.menu_padding-1.03, N),
            (1 - self.menu_padding, N))),) 
        
        glEnd()

        self._render_graph_info()

        for axis in range(3):
            time = 0

            glBegin(GL_LINE_STRIP)

            glColor3f(*self.colors[axis])
            
            for angvel in self.tbar.angvels:
                glVertex3f(M + time, angvel[axis] *\
                           self.graph_ratio *\
                           self.tbar.moment_inertia[axis]-0.5, 0)

                time += self.time_interval

            glEnd()

    def _render_graph_info1(self, color, m):
        glColor3f(*color)

        (*(glVertex3f(m, n, 0) for m,n in m),)

    """
    Screen._render_graph_info()
        renders the title and legend of the plot
    """
    def _render_graph_info(self):
        self._render_text((self.menu_padding-1 + 0.05, -0.1),\
                          "Angular momentum")
        
        glBegin(GL_LINES)

        N, Ñ = 0.85, 0.9

        (*(self._render_graph_info1(color, m) for color,m in
           ((YELLOW3F, ((N,-0.1), (Ñ,-0.1))),
            (GREEN3F, ((N,-0.17), (Ñ,-0.17))),
            (BLUE3F, ((N,-0.23), (Ñ,-0.23))))),)
        
        glEnd()

        (*(self._render_text((0.92, n), s_eje) for n, s_eje in\
           ((-0.11, "X"), (-0.18, "Y"), (-0.25, "Z"))),)

    """
    Screen._render_text(pos:tuple, text:str[, color:tuple])
        renders bitmap text
    """
    def _render_text(self, pos, text, color=(0, 0, 0, 255)):
        text_surface = self.font.render(text, True, color)
        text_data = tostring(text_surface, "RGBA", True)

        glRasterPos2d(*pos)
        glDrawPixels(text_surface.get_width(),\
                     text_surface.get_height(), GL_RGBA,\
                     GL_UNSIGNED_BYTE, text_data)

print('[*] Creating pygame display...')
init()
display, world, tbar = set_mode(WINDOW_SIZE, DOUBLEBUF | OPENGL),\
    World(), Tbar()

world.add_object(tbar)

screen = Screen(display=display, tbar=tbar)

print('[*] Initializing scene...')
world.init_scene()

frame, paused = 0, False
  
print('[*] Starting pygame loop...')
while frame < NUM_FRAMES:
    for event in get():
        if event.type == QUIT:
            print('[*] Quit event detected. Stopping pygame...')
            quit()
            exit(0)

        elif event.type == KEYUP:
            if event.key == K_SPACE:
                paused = not paused

                if paused: screen.show_paused_message()

    if not paused:
        world.clear()
        world.step()
        
        world.render()
        screen.render()

        frame += 1

        flip()
        wait(10)

print('[*] Simulation finished successfully')



