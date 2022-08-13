from manim import *
from colour import Color

g = 9.8
def wavelength(k):
    return 2*PI/k
def omega(k):
    return np.sqrt(2*PI*k*g)
def vph(k):
    return np.sqrt(2*PI/k*g)
def vgr(k):
    return np.sqrt(1/2*PI/k*g)


class V(Scene):
    config.background_color = Color("#244")

    def construct(self):
        title = Title("Angle of the wake")
        y0 = -2.5
        path = Line([10,y0,0],[-6,y0,0], stroke_color=RED)
        waves = []
        vphdots = []
        vgrdots = []
        for angle in np.arange(0,PI+.01,PI/10):
            waves += Line([6,y0,0],[6*np.cos(angle),y0+6*np.sin(angle),0], stroke_opacity=.5)
            vphdots += Dot([6*np.cos(angle),y0+6*np.sin(angle),0], color=GOLD)
            vgrdots += Dot([3+3*np.cos(angle),y0+3*np.sin(angle),0], color=YELLOW)
        vphdots = VGroup(*vphdots)
        vgrdots = VGroup(*vgrdots)

        angle = PI/3
        firstwave = Line([6,y0,0],[6*np.cos(angle),y0+6*np.sin(angle),0], stroke_opacity=.5)
        obsline = Line([-6,y0,0],[6*np.cos(angle),y0+6*np.sin(angle),0])

        circle_vph = Arc(radius = 6, arc_center = [0,y0,0], angle = PI, color=GOLD)
        circle_vgr = Arc(radius = 3, arc_center = [3,y0,0], angle = PI, color=YELLOW)

        tanline_long = Line([-6,y0,0],[8,y0+14/np.sqrt(8),0], stroke_color=WHITE)
        tanline = Line([-6,y0,0],[3-1,y0+np.sqrt(8),0], stroke_color=WHITE)
        trigline = Line([3,y0,0],[3-1,y0+np.sqrt(8),0], stroke_color=WHITE)
        traceline = Line([3,y0,0],[-6,y0,0], stroke_color=WHITE)

        b1 = Brace(traceline, direction=traceline.copy().rotate(PI / 2).get_unit_vector())
        b1text = b1.get_tex(r"3r")
        b1gr = VGroup(b1,b1text)

        b2 = Brace(trigline, direction=trigline.copy().rotate(-PI / 2).get_unit_vector())
        b2text = b2.get_tex(r"r")
        b2gr = VGroup(b2,b2text)

        a1 = RightAngle(trigline,tanline_long, quadrant=(-1, -1))
        a2 = Angle(traceline, tanline, quadrant=(-1, 1), radius = 1)
        a2label = MathTex(r"\theta")
        a2label.next_to(a2, RIGHT).shift(.1*UP)

        eqpos = 3*LEFT+UP
        eq1 = MathTex(r"\qquad \sin \theta =", r"\frac{r}{3r} = \frac13", r"\phantom{\approx 19.5^{\circ}}").move_to(eqpos)
        eq3 = MathTex(r"\theta = \arcsin", r"\frac{1}{3}", r"\phantom{\approx 19.5^{\circ}}").move_to(eqpos)
        eq4 = MathTex(r"\theta = \arcsin", r"\frac{1}{3}", r"\approx 19.5^{\circ}").move_to(eqpos)

        self.add(path)
        self.add(circle_vph)
        self.add(firstwave,obsline)
        self.wait(5)
        self.play(FadeOut(obsline,firstwave))
        self.play(Create(VGroup(*waves)), run_time=3)
        self.wait()
        self.play(Create(vphdots), run_time=2)
        self.wait(7)
        self.play( ReplacementTransform(vphdots,vgrdots), run_time=2 )
        self.wait(3)
        self.play(Create(circle_vgr))
        self.wait(5)

        self.play(Create(tanline_long))
        self.wait(3)
        self.play(FadeOut(*waves,*vphdots,*vgrdots,circle_vph))
        self.wait()
        self.add(tanline)
        self.play(Create(VGroup(trigline,traceline)),FadeIn(title))
        self.play(FadeOut(tanline_long))
        self.play(FadeIn(a1, a2, a2label, b1gr,b2gr))

        self.play(Write(eq1))
        self.wait()
        self.play( ReplacementTransform(eq1,eq3))
        self.play( ReplacementTransform(eq3,eq4))
        self.wait(5)

class Photo(Scene):
    def construct(self):
        l1 = Line([0,0,0],[-15,15*np.tan(np.arcsin(1/3)),0])
        l2 = Line([0,0,0],[-15,-15*np.tan(np.arcsin(1/3)),0])
        path = Line([-15,0,0],[0,0,0], stroke_color=RED)
        a1 = Angle(l1,path,quadrant=[1,-1], radius=2)
        a1label = MathTex(r"19.5^\circ")
        a1label.rotate(PI - 10*DEGREES).next_to(a1, LEFT).shift(.2*UP)
        vlines = VGroup(l1,path,l2,a1,a1label).rotate(PI).shift(9.5*RIGHT+.2*UP).set_opacity(.5)

        self.add(ImageMobject("elbe.png"))
        self.add(vlines)

class Photo2(Scene):
    def construct(self):
        l1 = Line([0,0,0],[15,15*np.tan(np.arcsin(1/3)),0])
        l2 = Line([0,0,0],[15,-15*np.tan(np.arcsin(1/3)),0])
        path = Line([15,0,0],[0,0,0], stroke_color=RED)
        a1 = Angle(l1,path,quadrant=[1,-1], radius=2, other_angle=True)
        a1label = MathTex(r"19.5^\circ")
        a1label.rotate(PI + 10*DEGREES).next_to(a1, RIGHT).shift(.2*UP)
        vlines = VGroup(l1,path,l2,a1,a1label).rotate(PI).shift(9.3*LEFT+.3*UP)#.set_opacity(.5)

        self.add(ImageMobject("rotterdam.png"))
        self.add(vlines)


