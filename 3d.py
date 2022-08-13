from manim import *
import scipy.integrate as integrate
import scipy.special as special
import scipy.signal as signal

#c = np.sqrt(9.5)/4 # approx sqrt(g) # previously .5 
# wave speed = c sqrt(labmda)
# period = c / sqrt(lambda)

# l = 2*PI/k
# v_ph = np.sqrt(2*PI/k*g) = np.sqrt(l*g)


g = 9.8/ 30 #gravitational constant in units of duck/second^2
# speed of duck is 1 
# time unit of ca .2s
# space unit of ca .5m
def wavelength(k):
    return 2*PI/k
def omega(k):
    return np.sqrt(2*PI*k*g)
def vph(k):
    return np.sqrt(2*PI/k*g)
def vgr(k):
    return np.sqrt(1/2*PI/k*g)

box = [[-15,2],[-7,7]]

time_step = .1
uv_cutoff = 1.0

def ellipsoid(A,B,C,x,y,z,color,res=8,**kwargs):
    func = lambda u,v:[-y - B*np.sin(u)*np.cos(v), x + A*np.cos(u)*np.cos(v), z + A*np.sin(v)]
    return Surface(func, u_range=[0,PI], v_range=[0,TAU], resolution = res, checkerboard_colors=[color], stroke_color=color, **kwargs)
duck = VGroup(
    ellipsoid(2,3,1.2,0,.75*PI,0, GRAY_BROWN, res = 16), #body
    ellipsoid(1,1.2,1,0,0,.5*PI, GREEN_E), #head
    ellipsoid(.4,1.3,.2,0,-.3*PI,.5*PI, YELLOW_B), #beak
    ellipsoid(.8,1.5,.2,0,1.4*PI,.1*PI, BLACK), #tail
    ellipsoid(.8,1.5,.2,0,-1.4*PI,.1*PI, BLACK, res=1, stroke_opacity=0, fill_opacity=0), #faketail
    ellipsoid(.15,.2,.15,.23*PI,-.15*PI,.65*PI, BLACK, res=4), #eye
    ellipsoid(.15,.2,.15,-.23*PI,-.15*PI,.65*PI, BLACK, res=4) #eye
).scale(.2).move_to([0,0,.5]).set_z_index(10)
  
#path = Surface(lambda u,v: [-20*u,.01*np.cos(v),.1+.01*np.sin(v)],u_range=[0,1], v_range=[0,TAU], stroke_color=RED, checkerboard_colors=[RED], stroke_width=2, resolution=16) #64
path = Line([box[0][0],0,0],[0,0,0], stroke_color=RED)

sea = Surface(lambda u,v: [u*box[0][0]+(1-u)*box[0][1],v*box[1][0]+(1-v)*box[1][1],0],u_range=[0,1], v_range=[0,1], stroke_color=BLUE, checkerboard_colors=[BLUE], fill_opacity=.5, stroke_opacity=.1, resolution=64)

l1 = Line([0,0,0],[-15,15*np.tan(np.arcsin(1/3)),0])
l2 = Line([0,0,0],[-15,-15*np.tan(np.arcsin(1/3)),0])
a1 = Angle(l1,path,quadrant=[1,-1], radius=3)
a1label = MathTex(r"19.5^\circ")
a1label.rotate(PI - 10*DEGREES).scale(1.8).next_to(a1, LEFT).shift(.2*UP)
vlines = VGroup(l1,path,l2,a1,a1label)

arrowkwds = {"max_stroke_width_to_length_ratio": 1e5, "max_tip_length_to_length_ratio": .7, "stroke_width": 6}
smallarrowkwds = {"max_stroke_width_to_length_ratio": 1e5, "max_tip_length_to_length_ratio": .3, "stroke_width": 3}

#def wave(start,end, c=2, phase=0, **kwargs):
#    diff = np.array(end)-np.array(start)
#    l = np.linalg.norm(diff)
#    normal = [0,0,.5]
#    func = lambda t: np.array(start) + t*diff + np.sin(t*l/wavelength[c] + phase)*normal
#    newphase = (l / wavelength[c]) % (2*np.pi)
#    return [ParametricFunction(func, t_range = np.array([0, 1]), **kwargs), newphase]


class Pebble(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=55 * DEGREES, theta=30 * DEGREES, frame_center=[0,0,-1])
        self.renderer.camera.light_source.move_to([-15,-7.5,15])

        t = ValueTracker(0)
        def updater():
            def surface(x,y):
                r = np.sqrt(x**2+y**2)*5
                a = .2
                tt = t.get_value()
                c = np.sqrt(2*a)*3
                func = lambda k: -np.cos(k*r - omega(k)*tt)*np.exp(-k**2/4/a)/c * 5/max(5,r)
                #func = lambda k: special.jv(0,k*r) * np.sin(-omega(k)*tt)
                z = integrate.quad(func, 0, np.inf, epsabs=1e-05, epsrel=1e-05)[0]
                return( [x,y, z] )
            waves = Surface(surface, v_range=[-6, 6], u_range=[-6, 6], resolution=64, checkerboard_colors=["#366"], stroke_color="#488") #128
            #pebble = Sphere(surface(0,0),radius=.1, resolution=4).set_color(RED)
            return waves #(VGroup(waves,pebble))

        waves = always_redraw(updater)
        self.add(waves)

        self.wait(2)
        self.play(t.animate.set_value(50),run_time=12,rate_func=linear)
        self.wait(2)

class PebbleMarked(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=55 * DEGREES, theta=30 * DEGREES, frame_center=[0,0,-1])
        self.renderer.camera.light_source.move_to([-15,-7.5,15])

        t = ValueTracker(0)
        def updater():
            tt = t.get_value()
            def radial(r):
                a = .2
                c = np.sqrt(2*a)*3
                func = lambda k: -np.cos(k*r - omega(k)*tt)*np.exp(-k**2/4/a)/c * 5/max(5,r)
                z = integrate.quad(func, 0, np.inf, epsabs=1e-05, epsrel=1e-05)[0]
                return z
            def surface(x,y):
                return( [x,y, radial(np.sqrt(x**2+y**2)*5)] )
               
            waves = Surface(surface, v_range=[-6, 6], u_range=[-6, 6], resolution=64, checkerboard_colors=["#366"], stroke_color="#488")#128
            if tt > 1:
                crests = signal.find_peaks([radial(r*5) for r in np.arange(0,6,.025)], height=.01)[0]
                circles = VGroup(*[Circle(radius=i*.025).move_to([0,0,radial(i*.025*5)]) for i in crests]).set_color(WHITE)
            else:
                circles = VGroup()    
            return(VGroup(waves,circles))

        waves = always_redraw(updater)
        self.add(waves)

        self.wait(2)
        self.play(t.animate.set_value(50),run_time=12,rate_func=linear) #50,12
        self.wait(2)


class Phase(ThreeDScene):
    def construct(self):
        #self.camera.frame_center=[-3., 0, 0.]
        self.set_camera_orientation(phi=60 * DEGREES, theta=30 * DEGREES, frame_center=[0,2,0])
        self.renderer.camera.light_source.move_to([-15,-7.5,15])

        self.add(sea,path)
        self.add_foreground_mobjects(duck)

        # set point of observation
        obs = Dot([-5,-3,0],fill_opacity = 1)
        obshighlight = Circle(radius=.2).move_to([-5,-3,0]).set_color(RED)

        self.wait(4)
        self.play(Create(obs))
        self.play(ShowPassingFlash(obshighlight, time_width=1))
        self.wait(2)

        duck.generate_target()

        for start in range(14,6,-1): #[8,6,4,2]: #starting positions for duck
            duck.target.shift(start*LEFT)
            if start == 14:
                self.play(MoveToTarget(duck))
            else:
                self.play(MoveToTarget(duck), run_time=0)

            # mark starting position / wave source
            source = Dot([-start,0,0],fill_opacity = 1, fill_color=BLUE)
            self.add(source)

            # get distance traveled
            diff = np.array(obs.get_center()) - np.array(source.get_center()) 
            dist = np.linalg.norm(diff)

            # compute wave number with correct phase velocity 
            # dist/start = v = np.sqrt(2*PI/k*g)
            vgr = dist/start
            vph = 2*vgr
            k = 2*PI*g/(vph)**2
            l = wavelength(k)
            normal = np.array([0,0,.5])

            # incrementally plot wave from source to obs
            def wave_updater():
                t = start+duck.get_center()[0] #time the duck has traveled (unit velocity)
                #func = lambda t: np.array([-start,0,0]) + t/start*diff + np.cos(t/start*dist/l)*normal
                func = lambda d: np.array([-start,0,0]) + d*diff/dist - np.sin(k*d - omega(k)*t)*normal
                wave = VGroup(
                    ParametricFunction(func, t_range = [t*vgr - l, t*vgr]),
                    Dot(np.array([-start,0,0]) + t*vgr*diff/dist - np.sin(k*t*vgr - omega(k)*t)*normal, fill_color=BLUE)
                )
                return(wave)
            wave = always_redraw(wave_updater)
            wave.set_shade_in_3d(True)
              
            # move duck to final time and draw wave
            self.add(wave)
            duck.target.shift(start*RIGHT)
            self.play(MoveToTarget(duck),run_time=start**2/14/3,rate_func=linear)

            # arrows to mark contribution of this particular wave
            arrow1 = Arrow(obs.get_center(),np.array(obs.get_center()) - np.sin(k*dist - omega(k)*start)*normal, **arrowkwds)
            arrow2 = Arrow(source.get_center(),np.array(source.get_center()) - np.sin(k*dist - omega(k)*start)*normal, **arrowkwds).set_color(RED)
            self.add(arrow1)
            self.wait(start/21)
            self.add(arrow2)
            self.wait(start/21)
            self.remove(arrow1)

            # remove wave before repeating with different source point
            self.remove(wave)

        self.wait(5)


class PhaseSample(ThreeDScene):
    def construct(self):
        #self.camera.frame_center=[-3., 0, 0.]
        self.set_camera_orientation(phi=60 * DEGREES, theta=30 * DEGREES, frame_center=[0,2,0])
        self.renderer.camera.light_source.move_to([-15,-7.5,15])

        self.add(sea,path)
        self.add_foreground_mobjects(duck)

        obs = Dot([0,-2,0],fill_opacity = .5)
        self.add(obs)

        def updater():
            arrows = []
            x1,x2,x3 = obs.get_center()
            tot = 0
            for t in np.arange(-20.0, 0.0, time_step):
                r = np.sqrt((x1-t)**2 + x2**2)
                vgr = r/t
                vph = 2*vgr
                k = 2*PI*g/(vph)**2
                phi = k*r - omega(k)*t
                #if l > uv_cutoff:
                #    influence = .05/r
                #else:
                influence = smooth(wavelength(k)/uv_cutoff)*.05/r
                if t>-15:
                    arrows += Arrow([t,0,0],[t, 0, 50*influence*np.sin(phi)],buff=0, color=RED, **smallarrowkwds)
                tot += influence*np.sin(phi)
            arrows += Arrow([x1,x2,0],[x1,x2,tot],buff=0, color=WHITE, **arrowkwds)
            #arrows += Tex(round(x1,2)).move_to([0,1,0])
            return(VGroup(*arrows))

        arrows = always_redraw(updater)
        self.add(arrows)
        self.wait()
        self.play(obs.animate.move_to([-5.8,-2,0]),run_time=7) #,rate_func=rush_from)
        self.wait()
        self.play(obs.animate.move_to([-6.2,-2,0]),run_time=2)
        self.wait()
        self.play(obs.animate.move_to([-6.6,-2,0]),run_time=2)
        self.wait(6)
        self.play(obs.animate.move_to([-15.0,-2,0]),run_time=3,rate_func=rush_into)
        self.wait(3)

class PhaseGrid(ThreeDScene):
    def construct(self):
        #self.camera.frame_center=[-3., 0, 0.]
        self.set_camera_orientation(phi=60 * DEGREES, theta=30 * DEGREES, frame_center=[0,2,0])
        self.renderer.camera.light_source.move_to([-15,-7.5,15])

        self.add(sea)
        self.add(duck)
        #self.begin_ambient_camera_rotation() #to avoid foreground/background isssues

        arrows = [VGroup()]
        for h in [1.0,.5,.25]:
            arrows += [[]]
            for x2 in np.arange(box[1][0]+h/2,box[1][1],h):
                for x1 in np.arange(box[0][0]+h/2,box[0][1],h):
                    tot = 0
                    for t in np.arange(-20.0, 0.0, time_step):
                        r = np.sqrt((x1-t)**2 + x2**2)
                        vgr = r/t
                        vph = 2*vgr
                        k = 2*PI*g/(vph)**2
                        phi = k*r - omega(k)*t
                        influence = .05/r
                        if wavelength(k) > uv_cutoff:
                            tot += influence*np.sin(phi)
                    #if not((x1+.5)**2+(x2+.5)**2 < 1):
                    arrows[-1] += [Arrow([x1,x2,0],[x1,x2,tot],buff=0, color=TEAL, **arrowkwds),Dot([x1,x2,0],fill_opacity=.1+abs(tot), color=WHITE)] 
            arrows[-1] = VGroup(*arrows[-1]).set_z_index(-1)
            #self.remove(arrows[-2])
            self.play(FadeOut(arrows[-2]),run_time=1)
            self.play(Create(arrows[-1], lag_ratio=0.5),run_time=4)
            self.wait()
        self.wait(5)

class GreenScreenDuck(ThreeDScene):
    def construct(self):
        self.camera.background_color='#00F'
        self.set_camera_orientation(phi=60 * DEGREES, theta=30 * DEGREES, frame_center=[0,2,0])
        self.renderer.camera.light_source.move_to([-15,-7.5,15])

        self.add(duck)


class Stationary(ThreeDScene):
    def construct(self):
        #self.camera.frame_center=[-3., 0, 0.]
        self.set_camera_orientation(phi=60 * DEGREES, theta=30 * DEGREES, frame_center=[0,2,0])
        self.renderer.camera.light_source.move_to([-15,-7.5,15])

        self.add(sea,path)
        self.add_foreground_mobjects(duck)
        duck.generate_target()

        # updater for a single stationary wave in direction of travel
        def wave_updater1():
            angle = 0
            vph = np.cos(angle)
            vgr = vph/2
            k = 2*PI*g/(vph)**2
            l = wavelength(k)
            normal = np.array([0,0,.5])

            waves = []
            start = 12
            t = start+duck.get_center()[0] #time the duck has traveled (unit velocity)
            func = lambda d: np.array([-start+d*np.cos(angle),d*np.sin(angle),0]) - np.sin(k*d - omega(k)*t)*normal
            waves += ParametricFunction(func, t_range = [0, t*vph], stroke_opacity=.5)
            return(VGroup(*waves))
        wave1 = always_redraw(wave_updater1)

        # updater for a single stationary wave at angle
        def wave_updater2():
            angle = -PI/6
            vph = np.cos(angle)
            vgr = vph/2
            k = 2*PI*g/(vph)**2
            l = wavelength(k)
            normal = np.array([0,0,.5])

            waves = []
            start = 12
            t = start+duck.get_center()[0] #time the duck has traveled (unit velocity)
            func = lambda d: np.array([-start+d*np.cos(angle),d*np.sin(angle),0]) - np.sin(k*d - omega(k)*t)*normal
            waves += ParametricFunction(func, t_range = [0, t*vph], stroke_opacity=.5)
            line = Line([duck.get_center()[0],0,0],[duck.get_center()[0]+12*np.sin(angle),-12*np.cos(angle),0])
            return(VGroup(*waves,line))
        wave2 = always_redraw(wave_updater2)

        # updater for wave front at angle
        def wave_updater3():
            angle = -PI/6
            vph = np.cos(angle)
            vgr = vph/2
            k = 2*PI*g/(vph)**2
            l = wavelength(k)
            normal = np.array([0,0,.5])

            waves = []
            for start in np.arange(0,12.1,1):
                t = start+duck.get_center()[0] #time the duck has traveled (unit velocity)
                if t > 0:
                    func = lambda d: np.array([-start+d*np.cos(angle),d*np.sin(angle),0]) - np.sin(k*d - omega(k)*t)*normal
                    waves += ParametricFunction(func, t_range = [0, t*vph], stroke_opacity=.5)
            line = Line([duck.get_center()[0],0,0],[duck.get_center()[0]+20*np.sin(angle),-20*np.cos(angle),0], stroke_color=WHITE)
            return(VGroup(*waves,line))
        wave3 = always_redraw(wave_updater3)

        # updater for wave front at different angle
        def wave_updater4():
            angle = -65*DEGREES
            vph = np.cos(angle)
            vgr = vph/2
            k = 2*PI*g/(vph)**2
            l = wavelength(k)
            normal = np.array([0,0,.5])

            waves = []
            for start in np.arange(0,12.1,1):
                t = start+duck.get_center()[0] #time the duck has traveled (unit velocity)
                if t > 0:
                    func = lambda d: np.array([-start+d*np.cos(angle),d*np.sin(angle),0]) - np.sin(k*d - omega(k)*t)*normal
                    waves += ParametricFunction(func, t_range = [0, t*vph], stroke_opacity=.5)
            line = Line([duck.get_center()[0],0,0],[duck.get_center()[0]+20*np.sin(angle),-20*np.cos(angle),0], stroke_color=WHITE)
            return(VGroup(*waves,line))
        wave4 = always_redraw(wave_updater4)

        # updater sweeping out stationary waves originating from single point
        source = Dot([-12,0,0],fill_opacity = 1, fill_color=BLUE)
        alpha = ValueTracker(0)
        def circle_updater():
            angle = alpha.get_value()
            vph = np.cos(angle)
            vgr = vph/2
            k = 2*PI*g/(vph)**2
            l = wavelength(k)
            normal = np.array([0,0,.5])

            start = 12
            t = 12 #time the duck has traveled (unit velocity)
            func = lambda d: np.array([-start+d*np.cos(angle),d*np.sin(angle),0]) - np.sin(k*d - omega(k)*t)*normal
            wave = ParametricFunction(func, t_range = [0, t*vph], stroke_opacity=.5)
            line = Line([duck.get_center()[0],0,0],[duck.get_center()[0]+20*np.sin(angle),-20*np.cos(angle),0], stroke_color=WHITE)
            return(VGroup(wave,line))
        circle = always_redraw(circle_updater)
        circle_arc = Arc(radius = 6, arc_center = [-6,0,0], angle = -PI, color=GOLD)

        self.wait(4)
        self.play(duck.animate.shift([-12,0,0]),run_time=1)
        self.wait(1)

        self.add(wave1)
        self.play(duck.animate.shift(12*RIGHT),run_time=8,rate_func=linear)
        self.wait(3)
        self.remove(wave1)
        duck.shift(12*LEFT)

        self.add(wave2)
        self.play(duck.animate.shift(12*RIGHT),run_time=17,rate_func=linear)
        self.wait(2)
        self.remove(wave2)
        duck.shift(12*LEFT)

        self.add(wave3)
        self.play(duck.animate.shift(12*RIGHT),run_time=6,rate_func=linear)
        self.wait()
        self.remove(wave3)
        duck.shift(12*LEFT)

        self.add(wave4)
        self.play(duck.animate.shift(12*RIGHT),run_time=6,rate_func=linear)
        self.wait(2)
        self.remove(wave4)

        self.add(circle, source)
        self.play(alpha.animate.set_value(-PI/2),run_time=4)
        self.add(circle_arc)
        self.play(alpha.animate.set_value(0),run_time=4)

        self.move_camera(phi=0, theta=0, gamma=-PI/2, frame_center=[-6,-.5,0], run_time=5, focal_distance=100) 
        self.play(alpha.animate.set_value(-PI/3),run_time=4)
        self.wait(3)


class Wake1(ThreeDScene):
    def construct(self):
        def height(x,y):
            tot = 0
            for t in np.arange(-20.0, 0.0, time_step):
                r = np.sqrt((x-t)**2 + y**2)
                vgr = r/t
                vph = 2*vgr
                k = 2*PI*g/(vph)**2
                phi = k*r - omega(k)*t
                influence = .01/r
                if wavelength(k) > uv_cutoff and r > .2:
                    tot += influence*np.sin(phi)
            return np.array([x,y,tot])
            
        wake = Surface(height, v_range=[-7, 7], u_range=[-15, 2], resolution=128,  checkerboard_colors=["#366"], stroke_color="#488")
        self.renderer.camera.light_source.move_to([-15,-7.5,15])
        self.set_camera_orientation(phi=60 * DEGREES, theta=30 * DEGREES, frame_center=[0,2,0])
        self.add(wake)
        self.add(duck)
        self.wait(5)
        self.move_camera(phi=65 * DEGREES, theta=195 * DEGREES, frame_center=[-8,2,-2], run_time=5)
        self.wait(10)


class Wake2(ThreeDScene):
    def construct(self):
       def height(x,y):
            tot = 0
            for t in np.arange(-20.0, 0.0, time_step):
                r = np.sqrt((x-t)**2 + y**2)
                vgr = r/t
                vph = 2*vgr
                k = 2*PI*g/(vph)**2
                phi = k*r - omega(k)*t
                influence = .01/r
                if wavelength(k) > uv_cutoff and r > .2:
                    tot += influence*np.sin(phi)
            return np.array([x,y,tot])

       wake = Surface(height, v_range=[-7, 7], u_range=[-15, 2], resolution=128, checkerboard_colors=["#366"], stroke_color="#488")
       self.renderer.camera.light_source.move_to([-15,-7.5,15])
       self.set_camera_orientation(phi=60 * DEGREES, theta=30 * DEGREES, frame_center=[0,2,0])
       self.add(wake)
       self.add(vlines)
       self.add_foreground_mobjects(duck)
       self.wait(10)
       self.remove(vlines)
       self.begin_3dillusion_camera_rotation(rate=.1)
       self.wait(30)
       self.stop_3dillusion_camera_rotation()
