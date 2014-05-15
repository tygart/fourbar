__version__ = '0.2.0'

from kivy.app import App
from kivy.clock import Clock
from kivy.uix.label import Label
from kivy.uix.floatlayout import FloatLayout
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.gridlayout import GridLayout
from kivy.uix.stacklayout import StackLayout
from kivy.uix.image import Image
from kivy.uix.button import Button
from kivy.uix.widget import Widget
from kivy.uix.slider import Slider
from kivy.properties import NumericProperty, ListProperty, StringProperty
from kivy.graphics import Color, Ellipse, Line, Rectangle
from kivy.metrics import dp, sp
from kivy.uix.popup import Popup
from kivy.uix.scrollview import ScrollView
from kivy.effects.scroll import ScrollEffect
from kivy.uix.dropdown import DropDown
from kivy.garden.navigationdrawer import NavigationDrawer
from kivy.garden.graph import Graph, MeshLinePlot, SmoothLinePlot
from kivy.utils import get_color_from_hex as rgb
import math
import json
from os.path import join, exists
from functools import partial


class FourBars(FloatLayout):

    type_bar = StringProperty('')
    mid = ListProperty([0, 0])

    def __init__(self, **kwargs):
        super(FourBars, self).__init__(**kwargs)
        self.pnt_a = Point(dp(0)+self.mid[0], dp(0)+self.mid[1], 'a')
        self.pnt_b = Point(dp(190)+self.mid[0], dp(0)+self.mid[1], 'b')
        self.pnt_d = Point(dp(30)+self.mid[0], dp(175)+self.mid[1], 'c')
        self.pnt_c = Point(dp(160)+self.mid[0], dp(175)+self.mid[1], 'd')
        self.pnts = [self.pnt_a, self.pnt_b, self.pnt_d, self.pnt_c]
        self.rad = dp(15.0)
        self.rocker_ind = None
        self.dir_test = None
        self.trace_c = []
        self.trace_d = []
        self.lengt = [1, 1, 1, 1]
        self.touch_test = True
        self.dir = 3
        self.hold_one = None
        self.data = {'map_x':[None,None,None,None],'map_y':[None,None,None,None],'name':[None,None,None,None]}
        self.omega_dc = []
        self.omega_ca = []
        self.omega_db = []
        self.pos_ca = []
        self.pos_dc = []
        self.pos_db = []
        self.max_time = 180

    def add_points(self):
        self.add_widget(self.pnt_a)
        self.add_widget(self.pnt_b)
        self.add_widget(self.pnt_c)
        self.add_widget(self.pnt_d)

    def graph_paper_draw(self, *args):

        self.graph_paper = Widget()

        delta = 35
        grid = self.build_grid(delta)
        grid_bold = self.build_grid(delta*5)

        with self.graph_paper.canvas:
            Color(0.8588, 0.8941, 0.60, 1)
            Rectangle(pos=self.pos, size=self.size)
            Color(0.3529, 0.3922, 0.102, 0.5)
            Line(points=grid, width=1)
            Line(points=grid_bold, width=dp(1.5))

        self.add_widget(self.graph_paper)

    def build_grid(self, delta):

        delta = dp(delta)
        max_grid = (dp(self.mid[0]*2), dp(self.mid[1]*2))

        grid = [0, 0]
        for x in range(0, int(max_grid[1]/(2*delta))+1, 1):
            grid.extend([0, grid[-1]+delta, max_grid[0], grid[-1]+delta, max_grid[0], grid[-1]+2*delta,
                         0, grid[-1]+2*delta])
        grid.extend([0, 0])
        for x in range(0, int(max_grid[0]/(2*delta))+1, 1):
            grid.extend([grid[-2]+delta, 0, grid[-2]+delta, max_grid[1],
                         grid[-2]+2*delta, max_grid[1], grid[-2]+2*delta, 0])

        return grid

    def reset(self, *args):

        min_x, min_y, max_x, max_y, hyp_pnt, hyp = None, None, None, None, None, None

        for pnt in self.pnts:
            test = (pnt.map_x**2+pnt.map_y**2)**0.5
            if hyp is None or test < hyp:
                hyp = test
                hyp_pnt = pnt
            if min_x is None or pnt.map_x <= min_x:
                min_x = pnt.map_x
            if max_x is None or pnt.map_x >= max_x:
                max_x = pnt.map_x
            if min_y is None or pnt.map_y <= min_y:
                min_y = pnt.map_y
            if max_y is None or pnt.map_y >= max_y:
                max_y = pnt.map_y

        dx = (max_x - min_x)/2.0
        dy = (max_y - min_y)/2.0
        zero_pnt = list([hyp_pnt.map_x, hyp_pnt.map_y])

        for pnt in self.pnts:
            pnt.map_x = (pnt.map_x + self.mid[0] - dx - zero_pnt[0])
            pnt.map_y = (pnt.map_y + self.mid[1] - dy - zero_pnt[1])

        self.lengt = self.calc_lengths()
        self.type_of(self.lengt)
        self.rocker_ind = None
        self.dir_test = None
        self.trace_c = []
        self.trace_d = []
        self.omega_dc = []
        self.omega_ca = []
        self.omega_db = []
        self.pos_ca = []
        self.pos_dc = []
        self.pos_db = []

    def switch_dir(self):
        self.dir *= -1

    def angle_domain(self, angle):
        if angle <= 0:
            angle += 360
        elif angle >= 360:
            angle -= 360
        return angle

    def trace_test(self, trace, x, y):
        if not trace:
            return True
        z = ((x - trace[-2])**2 + (y - trace[-1])**2)**0.5
        if z < dp(5.0):
            return False
        return True

    def update(self, dt):

        leng = self.lengt
        ang = self.calc_angles()
        old_d = (self.pnt_d.map_x, self.pnt_d.map_y, ang[1])
        old_c = (self.pnt_c.map_x, self.pnt_c.map_y, ang[0])

        if not self.type_bar == 'S + L > P + Q':
            self.graph_update(leng, ang)
        try:
            self.remove_widget(self.lines)
        except AttributeError as e:
            print(e)
            pass

        if self.type_bar == 'Drag Link' or self.type_bar == 'Crank-Rocker' or self.type_bar == 'Parallelogram Linkage':

            ang[0] += self.dir
            ang[0] = self.angle_domain(ang[0])

            self.pnt_c.map_x = leng[2]*math.cos(ang[0]*math.pi/180)+self.pnt_a.map_x
            self.pnt_c.map_y = leng[2]*math.sin(ang[0]*math.pi/180)+self.pnt_a.map_y

            alpha = leng[0]*math.cos(ang[2]*math.pi/180)-leng[2]*math.cos(ang[0]*math.pi/180)
            beta = leng[0]*math.sin(ang[2]*math.pi/180)-leng[2]*math.sin(ang[0]*math.pi/180)
            gamma = (leng[1]**2+alpha**2+beta**2-leng[3]**2)/(2*leng[1]*(alpha**2+beta**2)**0.5)

            if not self.type_bar == 'Parallelogram Linkage':
                ang[1] = self.fncn(alpha, beta, gamma, ang[1])
            else:
                ang[1] = ang[0]

            self.pnt_d.map_x = leng[1]*math.cos(ang[1]*math.pi/180)+self.pnt_b.map_x
            self.pnt_d.map_y = leng[1]*math.sin(ang[1]*math.pi/180)+self.pnt_b.map_y

            x_c = int(self.pnt_c.map_x+self.rad)
            y_c = int(self.pnt_c.map_y+self.rad)
            x_d = int(self.pnt_d.map_x+self.rad)
            y_d = int(self.pnt_d.map_y+self.rad)

            if len(self.trace_c) < 384 and self.trace_test(self.trace_c, x_c, y_c):
                self.trace_c.extend([x_c, y_c])
            if len(self.trace_d) < 384 and self.trace_test(self.trace_d, x_d, y_d):
                self.trace_d.extend([x_d, y_d])
                if len(self.trace_d) >= 4:
                    self.fast_trace(self.trace_d, old_d[2], ang[1], leng[1], self.pnt_b)

        elif self.type_bar == 'Rocker-Crank':

            ang[1] += self.dir
            ang[1] = self.angle_domain(ang[1])

            self.pnt_d.map_x = leng[1]*math.cos(ang[1]*math.pi/180)+self.pnt_b.map_x
            self.pnt_d.map_y = leng[1]*math.sin(ang[1]*math.pi/180)+self.pnt_b.map_y

            alpha = leng[0]*math.cos(ang[2]*math.pi/180)+leng[1]*math.cos(ang[1]*math.pi/180)
            beta = leng[0]*math.sin(ang[2]*math.pi/180)+leng[1]*math.sin(ang[1]*math.pi/180)
            gamma = (leng[3]**2-leng[2]**2-alpha**2-beta**2)/(2*leng[2]*(alpha**2+beta**2)**0.5)

            ang[0] = self.fncn(alpha, beta, gamma, ang[0])

            self.pnt_c.map_x = leng[2]*math.cos(ang[0]*math.pi/180)+self.pnt_a.map_x
            self.pnt_c.map_y = leng[2]*math.sin(ang[0]*math.pi/180)+self.pnt_a.map_y

            x_c = int(self.pnt_c.map_x+self.rad)
            y_c = int(self.pnt_c.map_y+self.rad)
            x_d = int(self.pnt_d.map_x+self.rad)
            y_d = int(self.pnt_d.map_y+self.rad)

            if len(self.trace_c) < 384 and self.trace_test(self.trace_c, x_c, y_c):
                self.trace_c.extend([x_c, y_c])
                if len(self.trace_c) >= 4:
                    self.fast_trace(self.trace_c, old_c[2], ang[0], leng[2], self.pnt_a)
            if len(self.trace_d) < 384 and self.trace_test(self.trace_d, x_d, y_d):
                self.trace_d.extend([x_d, y_d])

        elif self.type_bar == 'Double Rocker':

            if self.rocker_ind is None:
                self.rocker_ind = ang[3]

            self.rocker_ind += self.dir
            ang[3] = self.angle_domain(self.rocker_ind)

            alpha = leng[3]*math.cos(ang[3]*math.pi/180)-leng[0]*math.cos(ang[2]*math.pi/180)
            beta = leng[3]*math.sin(ang[3]*math.pi/180)-leng[0]*math.sin(ang[2]*math.pi/180)
            gamma = (leng[2]**2+alpha**2+beta**2-leng[1]**2)/(2*leng[2]*(alpha**2+beta**2)**0.5)

            ang[0] = self.fncn(alpha, beta, gamma, ang[0])

            self.pnt_c.map_x = leng[2]*math.cos(ang[0]*math.pi/180)+self.pnt_a.map_x
            self.pnt_c.map_y = leng[2]*math.sin(ang[0]*math.pi/180)+self.pnt_a.map_y

            self.pnt_d.map_x=leng[2]*math.cos(ang[0]*math.pi/180)+leng[3]*math.cos(ang[3]*math.pi/180)+self.pnt_a.map_x
            self.pnt_d.map_y=leng[2]*math.sin(ang[0]*math.pi/180)+leng[3]*math.sin(ang[3]*math.pi/180)+self.pnt_a.map_y
            ang = self.calc_angles()

            x_c = int(self.pnt_c.map_x+self.rad)
            y_c = int(self.pnt_c.map_y+self.rad)
            x_d = int(self.pnt_d.map_x+self.rad)
            y_d = int(self.pnt_d.map_y+self.rad)

            if len(self.trace_c) < 384 and self.trace_test(self.trace_c, x_c, y_c):
                self.trace_c.extend([x_c, y_c])
                if len(self.trace_c) >= 4:
                    self.fast_trace(self.trace_c, old_c[2], ang[0], leng[2], self.pnt_a)

            if len(self.trace_d) < 384 and self.trace_test(self.trace_d, x_d, y_d):
                self.trace_d.extend([x_d, y_d])
                if len(self.trace_d) >= 4:
                    self.fast_trace(self.trace_d, old_d[2], ang[1], leng[1], self.pnt_b)

        else:
            #self.type_bar == 'S + L > P + Q'
            pass

        self.add_lines()

    def fncn(self, alpha, beta, gamma, test):

        if alpha < 0:
            try:
                delta = math.atan(beta/alpha)*180/math.pi
            except ZeroDivisionError as e:
                print(e)
                delta = 0
        else:
            try:
                delta = math.atan(beta/alpha)*180/math.pi + 180
            except ZeroDivisionError as e:
                print(e)
                delta = 180
        try:
            theta1 = delta + math.acos(gamma)*180/math.pi
            theta2 = delta - math.acos(gamma)*180/math.pi
        except ValueError as e:
            print(e)
            print('tried to take arccos of %f' % gamma)
            print('other values', alpha, beta)
            return test

        theta1 = self.angle_domain(theta1)
        theta2 = self.angle_domain(theta2)

        if self.dir_test is None:
            if abs(theta1-test) < abs(theta2-test):
                self.dir_test = True
                return theta1
            else:
                self.dir_test = False
                return theta2
        elif self.dir_test is True:
            return theta1
        elif self.dir_test is False:
            return theta2

    def graph_update(self, leng, ang):

        omega = self.dir

        if self.type_bar == 'Drag Link' or self.type_bar == 'Crank-Rocker' or self.type_bar == 'Parallelogram Linkage':
            o_ca = omega
            o_dc = (-1*omega*leng[2]*math.sin((ang[1]-ang[0])*math.pi/180)) / \
                   (leng[3]*math.sin((ang[1]-ang[3])*math.pi/180))
            o_db = (omega*leng[2]*math.sin((ang[3]-ang[0])*math.pi/180)) / \
                   (leng[1]*math.sin((ang[3]-ang[1])*math.pi/180))

        elif self.type_bar == 'Rocker-Crank':
            o_ca = (omega*leng[1]*math.sin((ang[3]-ang[1])*math.pi/180)) / \
                   (leng[2]*math.sin((ang[3]-ang[0])*math.pi/180))
            o_dc = (-1*o_ca*leng[2]*math.sin((ang[1]-ang[0])*math.pi/180)) / \
                   (leng[3]*math.sin((ang[1]-ang[3])*math.pi/180))
            o_db = omega

        elif self.type_bar == 'Double Rocker':
            o_ca = (-1*omega*leng[3]*math.sin((ang[1]-ang[3])*math.pi/180)) / \
                   (leng[2]*math.sin((ang[1]-ang[0])*math.pi/180))
            o_dc = omega
            o_db = (o_ca*leng[2]*math.sin((ang[3]-ang[0])*math.pi/180)) / \
                   (leng[1]*math.sin((ang[3]-ang[1])*math.pi/180))

        else:
            #self.type_bar == 'S + L > P + Q':
            o_ca, o_dc, o_db = 0, 0, 0

        self.omega_ca.append(o_ca)
        self.omega_dc.append(o_dc)
        self.omega_db.append(o_db)
        self.pos_ca.append((leng[2], ang[0]))
        self.pos_dc.append((leng[3], ang[3]))
        self.pos_db.append((leng[1], ang[1]))

        if len(self.omega_ca) > self.max_time:
            self.omega_ca.pop(0)
            self.omega_dc.pop(0)
            self.omega_db.pop(0)
            self.pos_ca.pop(0)
            self.pos_dc.pop(0)
            self.pos_db.pop(0)

        if len(self.omega_ca) > self.max_time:
            Clock.schedule_once(self.reset, 0)

    def transition(self, dt):
        leng = self.calc_lengths()
        self.remove_widget(self.lines)
        self.type_of(leng)
        self.add_lines()

    def fast_trace(self, trace, old, new, leng, pnt):

        b_y = trace[-1]
        b_x = trace[-2]
        a_x = trace[-4]
        a_y = trace[-3]
        hyp = ((b_y - a_y)**2 + (b_x - a_x)**2)**0.5
        theta = math.acos((2*leng**2-hyp**2)/(2*leng**2))*180/math.pi
        c = 6
        a = c
        if (new - old < 0 or new - old > 180) and new - old > -180:
            a = -c

        alpha = old
        test = 0
        while True:
            alpha += a
            test += c
            if test >= theta:
                break
            alpha = self.angle_domain(alpha)
            trace.insert(-2, int(leng*math.cos(alpha*math.pi/180)+pnt.map_x+self.rad))
            trace.insert(-2, int(leng*math.sin(alpha*math.pi/180)+pnt.map_y+self.rad))

    def add_lines(self):
        a = int(self.rad)
        b = dp(2)
        bars = []
        for pnt in self.pnts:
            bars.append(pnt.map_x)
            bars.append(pnt.map_y)
        pnts = bars
        bars.append(self.pnt_a.map_x)
        bars.append(self.pnt_a.map_y)
        bars = [x+a for x in bars]

        self.lines = Widget()

        with self.lines.canvas:
            Color(0.2196, 0.6431, 0.8, 1)
            Line(points=self.trace_c, width=1)
            Color(0.8039, 0.2275, 0.3098, 1)
            Line(points=self.trace_d, width=1)
            Color(0.349, 0.349, 0.349, 1)
            Line(points=bars, width=b)

            Ellipse(pos=(pnts[0], pnts[1]), size=(a*2, a*2))
            Ellipse(pos=(pnts[2], pnts[3]), size=(a*2, a*2))
            Color(0.2, 0.2, 0.2, 1)
            Ellipse(pos=(pnts[0]+dp(4), pnts[1]+dp(4)), size=(a*2-dp(8), a*2-dp(8)))
            Ellipse(pos=(pnts[2]+dp(4), pnts[3]+dp(4)), size=(a*2-dp(8), a*2-dp(8)))
            Color(0.349, 0.349, 0.349, 1)
            Ellipse(pos=(pnts[6], pnts[7]), size=(a*2, a*2))
            Ellipse(pos=(pnts[4], pnts[5]), size=(a*2, a*2))
            Color(0.2196, 0.6431, 0.8, 1)
            Ellipse(pos=(pnts[6]+dp(4), pnts[7]+dp(4)), size=(a*2-dp(8), a*2-dp(8)))
            Color(0.8039, 0.2275, 0.3098, 1)
            Ellipse(pos=(pnts[4]+dp(4), pnts[5]+dp(4)), size=(a*2-dp(8), a*2-dp(8)))

        self.add_widget(self.lines)

    def calc_lengths(self):
        ba = ((self.pnt_b.get_pos()[0] - self.pnt_a.get_pos()[0])**2 +
              (self.pnt_b.get_pos()[1] - self.pnt_a.get_pos()[1])**2)**0.5
        db = ((self.pnt_d.get_pos()[0] - self.pnt_b.get_pos()[0])**2 +
              (self.pnt_d.get_pos()[1] - self.pnt_b.get_pos()[1])**2)**0.5
        ca = ((self.pnt_c.get_pos()[0] - self.pnt_a.get_pos()[0])**2 +
              (self.pnt_c.get_pos()[1] - self.pnt_a.get_pos()[1])**2)**0.5
        dc = ((self.pnt_d.get_pos()[0] - self.pnt_c.get_pos()[0])**2 +
              (self.pnt_d.get_pos()[1] - self.pnt_c.get_pos()[1])**2)**0.5

        return [ba, db, ca, dc]

    def calc_angles(self):

        theta_ca_zero = 180/math.pi*math.atan2((self.pnt_c.get_pos()[1]-self.pnt_a.get_pos()[1]),
                                               (self.pnt_c.get_pos()[0]-self.pnt_a.get_pos()[0]))
        if theta_ca_zero < 0:
            theta_ca_zero += 360
        theta_db_zero = 180/math.pi*math.atan2((self.pnt_d.get_pos()[1]-self.pnt_b.get_pos()[1]),
                                               (self.pnt_d.get_pos()[0]-self.pnt_b.get_pos()[0]))
        if theta_db_zero < 0:
            theta_db_zero += 360
        theta_ba_zero = 180/math.pi*math.atan2((self.pnt_b.get_pos()[1]-self.pnt_a.get_pos()[1]),
                                               (self.pnt_b.get_pos()[0]-self.pnt_a.get_pos()[0]))
        if theta_ba_zero < 0:
            theta_ba_zero += 360
        theta_dc_zero = 180/math.pi*math.atan2((self.pnt_d.get_pos()[1]-self.pnt_c.get_pos()[1]),
                                               (self.pnt_d.get_pos()[0]-self.pnt_c.get_pos()[0]))
        if theta_dc_zero < 0:
            theta_dc_zero += 360

        return [theta_ca_zero, theta_db_zero, theta_ba_zero, theta_dc_zero]

    def type_of(self, lengthy):

        lengths = list(lengthy)
        shortest = None
        for i in range(len(lengths)):
            if lengths[i] == min(lengths) or shortest is None:
                shortest = i

        ls = min(lengths)+max(lengths)
        lengths.remove(max(lengths))
        lengths.remove(min(lengths))
        pq = sum(lengths)

        if ls < pq and shortest == 0:
            self.type_bar = 'Drag Link'
        elif ls < pq and shortest == 1:
            self.type_bar = 'Rocker-Crank'
        elif ls < pq and shortest == 2:
            self.type_bar = 'Crank-Rocker'
        elif ls < pq and shortest == 3:
            self.type_bar = 'Double Rocker'
        elif ls == pq:
            self.type_bar = 'Parallelogram Linkage'
        else:
            self.type_bar = 'S + L > P + Q'

    def on_touch_move(self, touch):

        for pnt in self.pnts:
            if pnt.hold_me is True or (abs(pnt.map_x-touch.x) < 2*self.rad and abs(pnt.map_y-touch.y) < 2*self.rad):
                if self.hold_one is False or pnt.hold_me is True:
                    pnt.map_x, pnt.map_y = touch.x-self.rad, touch.y-self.rad
                    pnt.hold_me = True
                    self.hold_one = True

    def on_touch_down(self, touch):
        for pnt in self.pnts:
            pnt.hold_me = False
        self.hold_one = False
        self.touch_test = True

        if self.collide_point(*touch.pos):

            Clock.unschedule(self.update)
            self.rocker_ind = None
            self.dir_test = None
            self.trace_c = []
            self.trace_d = []
            self.touch_test = False
            self.omega_dc = []
            self.omega_ca = []
            self.omega_db = []
            self.pos_ca = []
            self.pos_dc = []
            self.pos_db = []

            Clock.schedule_interval(self.transition, 1.0/30.0)

    def on_touch_up(self, touch, *args):
        for pnt in self.pnts:
            pnt.hold_me = False
        self.hold_one = False

        if self.collide_point(*touch.pos) and self.touch_test is False:

            Clock.unschedule(self.transition)
            self.lengt = self.calc_lengths()

            self.touch_test = True
            Clock.schedule_interval(self.update, 1.0 / 30.0)


class Point(Widget):

    map_x, map_y = NumericProperty(0), NumericProperty(0)
    name = StringProperty(None)

    def __init__(self, x, y, name, **kwargs):
        super(Point, self).__init__(**kwargs)
        self.map_x = x
        self.map_y = y
        self.name = name
        self.hold_me = False

    def get_pos(self):
        return (self.map_x, self.map_y)


class RootWindow(FloatLayout):

    bar_leng = ListProperty([0, 0, 0, 0])
    bar_ls = ListProperty(['', '', '', ''])

    def __init__(self, **kwargs):
        super(RootWindow, self).__init__(**kwargs)
        self.save_path = None
        self.update_points = False
        self.plot_ca = MeshLinePlot(color=[0.2196, 0.6431, 0.8, 1])
        self.plot_db = MeshLinePlot(color=[0.8039, 0.2275, 0.3098, 1])
        self.plot_dc = MeshLinePlot(color=[0.7569, 0.4039, 0.7569, 1])
        self.y_max = 5
        self.y_min = -5
        self.y_tick = 1
        self.title = 'none'
        self.graph_type = 'none'

    def window_resize(self):
        self.ids.nav_drawer.toggle_state()

    def pause(self):

        Clock.unschedule(self.ids.fourbar.update)
        self.ids.fourbar.rocker_ind = None
        self.ids.fourbar.dir_test = None
        self.ids.fourbar.trace_c = []
        self.ids.fourbar.trace_d = []
        self.ids.fourbar.touch_test = True
        self.ids.fourbar.omega_dc = []
        self.ids.fourbar.omega_ca = []
        self.ids.fourbar.omega_db = []

        self.ids.fourbar.pos_ca = []
        self.ids.fourbar.pos_dc = []
        self.ids.fourbar.pos_db = []
        Clock.schedule_interval(self.ids.fourbar.transition, 1.0/30.0)

    def unpause(self):
        Clock.unschedule(self.ids.fourbar.transition)
        self.ids.fourbar.lengt = self.ids.fourbar.calc_lengths()
        Clock.schedule_interval(self.ids.fourbar.update, 1.0 / 30.0)

    def r_reset(self):
        self.legend_set_sandl()
        self.pause()
        Clock.schedule_once(self.ids.fourbar.reset, -1)
        self.unpause()

    def legend_set_sandl(self):

        self.bar_leng = self.ids.fourbar.lengt
        shortest = None
        longest = None
        for i in range(len(self.bar_leng)):
            if self.bar_leng[i] == min(self.bar_leng) or shortest is None:
                shortest = i
            if self.bar_leng[i] == max(self.bar_leng) or longest is None:
                longest = i

        for bar in range(len(self.bar_ls)):
            self.bar_ls[bar] = ''
        if not all(bar == self.bar_leng[0] for bar in self.bar_leng):
            self.bar_ls[shortest] = 'S'
            self.bar_ls[longest] = 'L'

    def r_speed_test(self):

        if abs(self.ids.fourbar.dir) == 1 and self.ids.fourbar.max_time != 390:
            self.ids.fourbar.max_time = 390
            self.update_points = False
            self.graph_kill()
            Clock.schedule_once(partial(self.graph_build, self.graph_type), -1)
        elif abs(self.ids.fourbar.dir) == 2 and self.ids.fourbar.max_time != 270:
            self.ids.fourbar.max_time = 270
            self.update_points = False
            self.graph_kill()
            Clock.schedule_once(partial(self.graph_build, self.graph_type), -1)
        elif abs(self.ids.fourbar.dir) >= 3 and self.ids.fourbar.max_time != 180:
            self.ids.fourbar.max_time = 180
            self.update_points = False
            self.graph_kill()
            Clock.schedule_once(partial(self.graph_build, self.graph_type), -1)

    def update_leng(self, *args):

        self.legend_set_sandl()
        self.r_speed_test()

        if self.update_points is True and len(self.ids.fourbar.omega_ca) > 90 and self.graph_type is not 'none':

            if self.graph_type is 'velocity':

                self.test_y(self.calc_y_axis())
                self.plot_ca.points = [(x, y) for x, y in enumerate(self.ids.fourbar.omega_ca)]
                self.plot_db.points = [(x, y) for x, y in enumerate(self.ids.fourbar.omega_db)]
                self.plot_dc.points = [(x, y) for x, y in enumerate(self.ids.fourbar.omega_dc)]

            elif self.graph_type is 'delta_x':

                self.test_y_pos()
                self.plot_ca.points=[(x, y[0]*math.cos(y[1]*math.pi/180)) for x,y in enumerate(self.ids.fourbar.pos_ca)]
                self.plot_db.points=[(x, y[0]*math.cos(y[1]*math.pi/180)) for x,y in enumerate(self.ids.fourbar.pos_db)]
                self.plot_dc.points=[(x, y[0]*math.cos(y[1]*math.pi/180)) for x,y in enumerate(self.ids.fourbar.pos_dc)]

            elif self.graph_type is 'delta_y':

                self.test_y_pos()
                self.plot_ca.points=[(x, y[0]*math.sin(y[1]*math.pi/180)) for x,y in enumerate(self.ids.fourbar.pos_ca)]
                self.plot_db.points=[(x, y[0]*math.sin(y[1]*math.pi/180)) for x,y in enumerate(self.ids.fourbar.pos_db)]
                self.plot_dc.points=[(x, y[0]*math.sin(y[1]*math.pi/180)) for x,y in enumerate(self.ids.fourbar.pos_dc)]

    def test_y(self, y):

        if y[0] < self.y_min or y[1] > self.y_max or abs(y[0] - self.y_min) > 10 or abs(y[1] - self.y_max) > 10 \
                or self.y_tick == 0:
            self.update_points = False
            self.graph_kill()
            Clock.schedule_once(partial(self.graph_build, self.graph_type), 1)

    def test_y_pos(self):

        if self.y_max != int(max(self.ids.fourbar.lengt)):
            self.update_points = False
            self.graph_kill()
            Clock.schedule_once(partial(self.graph_build, self.graph_type), 1)

    def resize(self):

        try:
            self.ids.fourbar.remove_widget(self.ids.fourbar.graph_paper)
            Clock.schedule_once(self.ids.fourbar.graph_paper_draw, 0)
            self.graph_kill()
            self.update_points = False
            if self.graph_type is not 'none':
                Clock.schedule_once(partial(self.graph_build, self.graph_type), 0)
            self.r_reset()
            self.ids.nav_drawer.adjust_sidepanel()
        except AttributeError as e:
            print(e)
            pass

    def switch(self):
        self.pause()
        self.ids.fourbar.switch_dir()
        self.unpause()

    def start_step(self, dire):
        if dire < 0 and self.ids.fourbar.dir > 0:
            self.ids.fourbar.switch_dir()
        elif dire > 0 and self.ids.fourbar.dir < 0:
            self.ids.fourbar.switch_dir()
        self.pause()
        self.unpause()

    def load_points(self):
        if not exists(self.points_fn):
            return
        with open(self.points_fn, 'rb') as fd:
            data = json.load(fd)
        self.ids.fourbar.data = data

        for i, pnt in enumerate(self.ids.fourbar.pnts):
            pnt.name = self.ids.fourbar.data['name'][i]
            pnt.map_x = self.ids.fourbar.data['map_x'][i]
            pnt.map_y = self.ids.fourbar.data['map_y'][i]
        self.r_reset()

    def save_points(self):
        for i,pnt in enumerate(self.ids.fourbar.pnts):
            self.ids.fourbar.data['name'][i] = pnt.name
            self.ids.fourbar.data['map_x'][i] = pnt.map_x
            self.ids.fourbar.data['map_y'][i] = pnt.map_y

        with open(self.points_fn, 'wb') as fd:
            json.dump(self.ids.fourbar.data, fd)

    @property
    def points_fn(self):
        return join(self.save_path, 'fourbar.json')

    def example(self, beta, condit=False, *args):

        pnts = self.ids.fourbar.pnts  # [a,b,d,c]

        if beta is 'dc':
            pnts[0].map_x, pnts[0].map_y = dp(0), dp(0)
            pnts[1].map_x, pnts[1].map_y = dp(75), dp(0)
            pnts[3].map_x, pnts[3].map_y = dp(0), dp(90)
            pnts[2].map_x, pnts[2].map_y = dp(170), dp(175)
        elif beta is 'dr':
            pnts[0].map_x, pnts[0].map_y = dp(0), dp(0)
            pnts[1].map_x, pnts[1].map_y = dp(190), dp(0)
            pnts[3].map_x, pnts[3].map_y = dp(30), dp(175)
            pnts[2].map_x, pnts[2].map_y = dp(160), dp(175)
        elif beta is 'cr':
            pnts[0].map_x, pnts[0].map_y = dp(0), dp(0)
            pnts[1].map_x, pnts[1].map_y = dp(200), dp(10)
            pnts[3].map_x, pnts[3].map_y = dp(0), dp(125)
            pnts[2].map_x, pnts[2].map_y = dp(190), dp(200)
        elif beta is 'rc':
            pnts[0].map_x, pnts[0].map_y = dp(20), dp(200)
            pnts[1].map_x, pnts[1].map_y = dp(175), dp(185)
            pnts[3].map_x, pnts[3].map_y = dp(0), dp(60)
            pnts[2].map_x, pnts[2].map_y = dp(160), dp(65)
        elif beta is 'p':
            pnts[0].map_x, pnts[0].map_y = dp(0), dp(0)
            pnts[1].map_x, pnts[1].map_y = dp(100), dp(0)
            pnts[3].map_x, pnts[3].map_y = dp(0), dp(150)
            pnts[2].map_x, pnts[2].map_y = dp(100), dp(150)
        elif beta is 'reset':
            pnts[0].map_x, pnts[0].map_y = dp(0), dp(0)
            pnts[1].map_x, pnts[1].map_y = dp(100), dp(0)
            pnts[3].map_x, pnts[3].map_y = dp(0), dp(100)
            pnts[2].map_x, pnts[2].map_y = dp(100), dp(100)

        if condit is True:

            self.ids.nav_drawer.adjust_sidepanel()
            self.ids.fourbar_scroll.scroll_y = 0.5
            self.ids.fourbar_scroll.scroll_x = 0.5
            Clock.schedule_once(self.intro, 0)
            Clock.schedule_interval(self.update_leng, 1.0/30.0)

        self.r_reset()

    def intro(self, *args):
        p = GrashofPopup()
        p.open()
        self.ids.menu_init.bind(on_release=self.ids.menu_graph.open)

    def round_up(self, x):
        return int(math.ceil(x / 5.0)) * 5

    def round_down(self, x):
        return int(math.floor(x / 5.0)) * 5

    def graph_build(self, alpha, *args):

        self.graph_kill()
        self.graph_type = alpha

        if alpha is 'velocity':
            self.update_points = True
            x = self.calc_y_axis()
            if x[2] == 0:
                self.y_tick = 0
                return

            self.y_min, self.y_max, self.y_tick, self.title = x[0], x[1], x[2], 'Angular Velocity'
            self.graphy()
            self.plot_ca.points = [(x, y) for x, y in enumerate(self.ids.fourbar.omega_ca)]
            self.plot_db.points = [(x, y) for x, y in enumerate(self.ids.fourbar.omega_db)]
            self.plot_dc.points = [(x, y) for x, y in enumerate(self.ids.fourbar.omega_dc)]
            self.graph.add_plot(self.plot_ca)
            self.graph.add_plot(self.plot_db)
            self.graph.add_plot(self.plot_dc)
            self.add_widget(self.graph)

        if alpha is 'delta_x':
            self.update_points = True

            self.y_max = int(max(self.ids.fourbar.lengt))
            self.y_min = -self.y_max
            self.y_tick = int(2.0*self.y_max/5.0)
            self.title = 'Horizontal Translation'

            self.graphy()
            self.plot_ca.points = [(x, y[0]*math.cos(y[1]*math.pi/180)) for x, y in enumerate(self.ids.fourbar.pos_ca)]
            self.plot_db.points = [(x, y[0]*math.cos(y[1]*math.pi/180)) for x, y in enumerate(self.ids.fourbar.pos_db)]
            self.plot_dc.points = [(x, y[0]*math.cos(y[1]*math.pi/180)) for x, y in enumerate(self.ids.fourbar.pos_dc)]
            self.graph.add_plot(self.plot_ca)
            self.graph.add_plot(self.plot_db)
            self.graph.add_plot(self.plot_dc)
            self.add_widget(self.graph)

        if alpha is 'delta_y':
            self.update_points = True

            self.y_max = int(max(self.ids.fourbar.lengt))
            self.y_min = -self.y_max
            self.y_tick = int(2.0*self.y_max/5.0)
            self.title = 'Vertical Translation'

            self.graphy()
            self.plot_ca.points = [(x, y[0]*math.sin(y[1]*math.pi/180)) for x,y in enumerate(self.ids.fourbar.pos_ca)]
            self.plot_db.points = [(x, y[0]*math.sin(y[1]*math.pi/180)) for x,y in enumerate(self.ids.fourbar.pos_db)]
            self.plot_dc.points = [(x, y[0]*math.sin(y[1]*math.pi/180)) for x,y in enumerate(self.ids.fourbar.pos_dc)]
            self.graph.add_plot(self.plot_ca)
            self.graph.add_plot(self.plot_db)
            self.graph.add_plot(self.plot_dc)
            self.add_widget(self.graph)

    def calc_y_axis(self):

        if len(self.ids.fourbar.omega_ca) < 90 or not self.ids.fourbar.omega_ca:
            return 0, 0, 0

        data = [self.ids.fourbar.omega_ca, self.ids.fourbar.omega_db, self.ids.fourbar.omega_dc]
        y_max = 0
        y_min = 0
        for dat in data:

            yma = max(dat)
            ymi = min(dat)
            if ymi < y_min:
                y_min = ymi
            if yma > y_max:
                y_max = yma

        y_max = self.round_up(y_max)
        y_min = self.round_down(y_min)
        delta = abs(y_max - y_min)
        y_tick = delta/5

        return y_min, y_max, y_tick

    def graphy(self):

        graph_theme = {'label_options': {'color': rgb('595959'), 'bold': False}, 'background_color': rgb('DBE49A'),
                       'tick_color': rgb('999999'), 'border_color': rgb('808080')}

        self.graph = CustomGraph(pos=(self.width-dp(410), self.height - dp(347)), size_hint = (None, None),
                                 size = (dp(410), dp(300)), xlabel='Time', ylabel=self.title, x_ticks_minor=5,
                                 x_ticks_major=30, y_ticks_major=self.y_tick, y_grid_label=True, x_grid_label=True,
                                 padding=10, xlog=False, ylog=False, x_grid=True, y_grid=True, xmin=0,
                                 xmax=self.ids.fourbar.max_time, ymin=self.y_min, ymax=self.y_max, draw_border = True,
                                 **graph_theme)

    def graph_kill(self):

        try:
            self.graph.remove_plot(self.plot_ca)
            self.remove_widget(self.graph)
        except AttributeError as e:
            print(e)
            pass

    def speed(self):

        if self.ids.fourbar.dir > 0:
            self.ids.fourbar.dir = self.ids.speed.value
        else:
            self.ids.fourbar.dir = -1*self.ids.speed.value
        self.pause()
        self.unpause()

class CustomGraph(Graph):

    def __init__(self, **kwargs):
        super(CustomGraph, self).__init__(**kwargs)

        with self._fbo:
            Color(0.349, 0.349, 0.349, 1)
            Rectangle(pos=(dp(0), dp(0)), size=(dp(5), dp(300)))
            Rectangle(pos=(dp(0), dp(295)), size=(dp(410), dp(5)))
            Rectangle(pos=(dp(405), dp(295)), size=(dp(5), -dp(295)))
            Rectangle(pos=(dp(0), dp(0)), size=(dp(410), dp(5)))


class SideWindow(BoxLayout):

    def __init__(self, **kwargs):
        super(SideWindow, self).__init__(**kwargs)
        self.slider = 1

class BaseBar(Widget):

    def __init__(self, **kwargs):
        super(BaseBar, self).__init__(**kwargs)


class CustomSlider(Slider):

    def __init__(self, **kwargs):
        super(CustomSlider, self).__init__(**kwargs)


class DummyEffect(ScrollEffect):
    pass


class DummyScroll(ScrollView):

    def __init__(self, **kwargs):
        super(DummyScroll, self).__init__(**kwargs)
        self.effect_cls = DummyEffect


class DummyNavDrawer(NavigationDrawer):

    def __init__(self, **kwargs):
        super(DummyNavDrawer, self).__init__(**kwargs)

    def adjust_sidepanel(self):
        self.side_panel_width = dp(200)
        self.separator_image = 'data/navigationdrawer_gradient_ltor.png'


class CustomButton(Button):

    def __init__(self, **kwargs):
        super(CustomButton, self).__init__(**kwargs)


class GrashofPopup(Popup):

    grashof_con_text = StringProperty('')

    def __init__(self, **kwargs):
        super(GrashofPopup, self).__init__(**kwargs)
        self.size_hint_y = 0.9
        self.size_hint_x = 0.9
        self.grashof_con_text = u'If the sum of the shortest and longest link of a planar quadrilateral linkage ' + \
                                u'is less than or equal to the sum of the remaining two links, then the shortest ' + \
                                u'link can rotate fully with respect to a neighboring link. In other words, the ' + \
                                u'condition is satisfied if S+L ' + unichr(8804) + u' P+Q where S is the shortest ' + \
                                u'link, L is the longest, and P and Q are the other links.\n\nThis app allows you ' +\
                                u'to examine the motion of any four-bar linkage that satisfies Grashof\'s ' +\
                                u'condition.  You can adjust the lengths of the bars by moving the end points.\n'

    def on_size(self, instance, value):

        if self.height < self.width:
            self.ids.cases_pic.source = 'data/Linkage_four_bar.png'
        else:
            self.ids.cases_pic.source = 'data/Linkage_four_bar_square.png'


class FourBarApp(App):

    def build(self):

        root = RootWindow()
        root.save_path = self.user_data_dir
        root.ids.fourbar.add_points()
        root.ids.fourbar.add_lines()
        root.ids.nav_drawer.toggle_state()
        Clock.schedule_once(partial(root.example, 'dr', True), 0)
        Clock.schedule_once(root.ids.fourbar.graph_paper_draw, -1)

        return root


if __name__ == '__main__':
    FourBarApp().run()
