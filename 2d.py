from manim import *

g = 9.8
def wavelength(k):
    return 2*PI/k
def omega(k):
    return np.sqrt(2*PI*k*g)
def vph(k):
    return np.sqrt(2*PI/k*g)
def vgr(k):
    return np.sqrt(1/2*PI/k*g)

# Method to write several (aligned) lines. In hindsight not a great idea, in case you later want to manipulate/highlight individual lines
def write_lines(*content, align_left=True, align_right=False, indent=0, linespace = .7):
    out = [Tex(content[0])]
    for i in range(len(content))[1:]:
        if align_left and not(content[i]==""):
            if indent == 0:
                out += [Tex(content[i]).shift(linespace*i*DOWN).align_to(out[0],LEFT)]
            else:
                out += [Tex(content[i]).shift(linespace*i*DOWN).align_to(out[0],LEFT).shift(indent[i]*RIGHT)]
        elif align_right and not(content[i]==""):
            if indent == 0:
                out += [Tex(content[i]).shift(linespace*i*DOWN).align_to(out[0],RIGHT)]
            else:
                out += [Tex(content[i]).shift(linespace*i*DOWN).align_to(out[0],RIGHT).shift(indent[i]*LEFT)]
        else:
            out += [Tex(content[i]).shift(linespace*i*DOWN)]
    return VGroup(*out)

class SineWave(Scene):
    def construct(self):
        title = Title(r"Speed of wave depends on wavelength $\lambda$")

        ax1 = Axes(
            x_range=[-5, 15], y_range=[-1, 1], axis_config={"include_tip": False},
            x_length=15, y_length=.7, y_axis_config = {"stroke_width": 0}#x_length=6,
        )
        ax2 = ax1.copy()
        ax3 = ax1.copy()
      
        t = ValueTracker(0)

        k1 = 4.0
        k2 = 2.0
        k3 = 1.0

        def get_plot1():
            l = Line(ax1.coords_to_point(-3/2*PI/k1 + omega(k1)/k1*t.get_value(),1),ax1.coords_to_point(1/2*PI/k1 + omega(k1)/k1*t.get_value(),1), stroke_width=0)
            b = Brace(l, direction=[0,1,0], buff=0)
            bt = b.get_tex(r"\lambda")
            out = VGroup(
                ax1.plot(lambda x: np.sin(k1*x - omega(k1)*t.get_value()), color=BLUE),
                l,b,bt
            )
            return out 
        graph1 = always_redraw( get_plot1 )

        def get_plot2():
            l = Line(ax2.coords_to_point(-3/2*PI/k2 + omega(k2)/k2*t.get_value(),1),ax2.coords_to_point(1/2*PI/k2 + omega(k2)/k2*t.get_value(),1), stroke_width=0)
            b = Brace(l, direction=[0,1,0], buff=0)
            bt = b.get_tex(r"\lambda")
            out = VGroup(
                ax2.plot(lambda x: np.sin(k2*x - omega(k2)*t.get_value()), color=BLUE),
                l,b,bt
            )
            return out 
        graph2 = always_redraw( get_plot2 )

        def get_plot3():
            l = Line(ax3.coords_to_point(-3/2*PI/k3 + omega(k3)/k3*t.get_value(),1),ax3.coords_to_point(1/2*PI/k3 + omega(k3)/k3*t.get_value(),1), stroke_width=0)
            b = Brace(l, direction=[0,1,0], buff=0)
            bt = b.get_tex(r"\lambda")
            out = VGroup(
                ax3.plot(lambda x: np.sin(k3*x - omega(k3)*t.get_value()), color=BLUE),
                l,b,bt
            )
            return out 
        graph3 = always_redraw( get_plot3 )

        plot1 = VGroup(ax1, graph1)
        plot2 = VGroup(ax2, graph2)
        plot3 = VGroup(ax3, graph3)
        plots = VGroup(plot1,plot2,plot3).arrange(DOWN).shift(DOWN/2)

        
        self.add(title)
        self.add(plots)
        self.play( t.animate.set_value(2.5) , run_time=15, rate_func=linear )
        self.wait()

class SineWave2(Scene):
    def construct(self):
        title = Title(r"Formula for a travelling sine wave")

        ax1 = Axes(
            x_range=[0, 16, PI/2], y_range=[-2, 2], 
            x_length=10, y_length=3, #x_length=6, axis_config={"include_tip": False},
        ).shift(1.5*RIGHT+UP)
        axlabel = ax1.get_axis_labels(x_label="x", y_label="")

        k = ValueTracker(1)
        k1 = k.get_value()
        k2 = 3
        t = ValueTracker(0)

        def get_plot1():
            k1 = k.get_value()
            #l = Line(ax1.coords_to_point(1/2*PI/k1 + omega(k1)/k1*t.get_value(),1),ax1.coords_to_point(5/2*PI/k1 + omega(k1)/k1*t.get_value(),1), stroke_width=0)
            #b = Brace(l, direction=[0,1,0], buff=0)
            #bt = b.get_tex(r"\lambda")
            out = VGroup(
                ax1.plot(lambda x: np.sin(k1*x - omega(k1)*t.get_value()), color=BLUE),
                #l,b,bt
            )
            return out 
        graph1 = always_redraw( get_plot1 )
        plot1 = VGroup(ax1, graph1)

        dots = []
        labels = [r"kx = 90^\circ", r"kx = 90^\circ + 360^\circ", r"kx = 90^\circ + 720^\circ"]
        for i in range(3):
            dot = Dot(ax1.coords_to_point(1/2*PI/k1 + i*2*PI/k1,1))
            label = MathTex(labels[i]).next_to(dot,UP,.2).scale(.8)
            dots += [VGroup(dot,label)]
        dotsgroup = VGroup(*dots)

        text1 = MathTex(r"\sin",r"(k x)").shift(5.5*LEFT+1.5*UP).set_color(BLUE)
        text1b = MathTex(r"\sin",r"(k x - a)").shift(5.5*LEFT+1.5*UP).set_color(BLUE)
        text1c = MathTex(r"\sin",r"(k x - \omega t)").shift(5.5*LEFT+1.5*UP).set_color(BLUE)
        framebox1 = SurroundingRectangle(text1c, buff = .25).set_color(WHITE)
        highlight1 = Circle(radius=.3).set_color(RED).move_to(text1[1].get_center()).shift(.2*LEFT)
        highlight1b = Circle(radius=.4).set_color(RED).move_to(text1[1].get_center()).shift(.4*RIGHT)

        text2 = write_lines(
            r"Wave repeats itself when $kx$ increases by $360^\circ$",
            r"$x$ increases by $\displaystyle\frac{360^\circ}{k}$",
            linespace = 1,
            align_left = False,
            align_right = True
        ).shift(1*DOWN + 1.5*LEFT)
        text2b = MathTex(r"\Rightarrow").next_to(text2[1],RIGHT)
        text2c = MathTex(r"\lambda = \frac{360^\circ}{k}").next_to(text2b,RIGHT).set_color(YELLOW)
            
        def get_text3():
            k1 = k.get_value()
            return MathTex(r"k = " + str('{:.1f}'.format(k1)) ).shift(5.5*LEFT+.5*UP).set_color(YELLOW)
        text3 = always_redraw( get_text3 )
        def get_text3b():
            a = t.get_value()*omega(k.get_value())
            return MathTex(r"a = " + str('{:.1f}'.format(a)) ).shift(5.5*LEFT+.5*UP).set_color(YELLOW)
        text3b = always_redraw( get_text3b )
        def get_text3c():
            tt = t.get_value()
            out = VGroup(
                MathTex(r"t = " + str('{:.1f}'.format(tt)) ).shift(5.5*LEFT+.5*UP),
                #MathTex(r"\omega = " + str('{:.1f}'.format(omega(k2))) ).shift(5.5*LEFT+.2*DOWN)
            ).set_color(YELLOW)
            return out
        text3c = always_redraw( get_text3c )

        text4 = MathTex(r"kx - \omega t", r"= \mathrm{const}").move_to(RIGHT+2*DOWN).set_color(YELLOW)
        text4a = MathTex(r"kx - \omega t", r"= 0").move_to(RIGHT+2*DOWN).set_color(YELLOW)
        text4b = MathTex(r"\frac{x}{t}" , r" = \frac{\omega}{k}").move_to(RIGHT+2*DOWN).set_color(YELLOW)
        text4c = MathTex(r"v", r" = \frac{\omega}{k}").move_to(RIGHT+2*DOWN).set_color(YELLOW)
        
        text5 = MathTex(r"v ").move_to(4*LEFT+2*DOWN).set_color(YELLOW)
        text5a = MathTex(r"\frac{\omega}{k}", r" = \sqrt{\lambda g}").move_to(4*LEFT+2*DOWN).set_color(YELLOW)
        text5b = MathTex(r"\frac{\omega}{k}", r" = \sqrt{\frac{360^\circ \cdot g}{k}}").move_to(4*LEFT+2*DOWN).set_color(YELLOW)
        text5c = MathTex(r"\omega = \sqrt{360^\circ \cdot g \cdot k").move_to(4*LEFT+2*DOWN).set_color(YELLOW)

        self.play( Write(title) )
        self.wait()
        self.add( ax1,axlabel )
        self.play( Create(graph1), run_time=2, rate_func=linear )
        self.play( Write(text1) )
        self.wait(2)
        self.play( ShowPassingFlash(highlight1,time_width=1) )
        self.wait(3)
        self.play( Write(dotsgroup), run_time=3 )
        self.wait(1)

        # determine wavelength
        self.play( Write(text2[0]) )
        self.wait(4)
        self.play( Write(text2[1]) )
        self.wait(6)
        self.play( Write(VGroup(text2b, text2c)) )
        self.play( FadeOut(dotsgroup) )
        self.wait(1)
        self.play( FadeIn(text3) )
        self.play( k.animate.set_value(k2) , run_time=4 )
        self.remove(text2[0],text2[1],text2b)
        self.remove(text3)
        self.wait(2)

        # add phase shift
        self.play( ReplacementTransform(text1,text1b) )
        self.play( ShowPassingFlash(highlight1b,time_width=1) )
        self.wait(2)
        self.play( FadeIn(text3b) )
        self.play( t.animate.set_value(1) , run_time=4, rate_func=linear )
        self.remove(text3b)

        # add time dependence
        self.play( ReplacementTransform(text1b,text1c) )
        self.wait(3)
        self.play( Create(framebox1) )
        self.play( FadeIn(text3c) )
        self.play( t.animate.set_value(5) , run_time=10, rate_func=linear )
        self.wait(2)
        self.remove(text3c)

        # determine speed
        self.play( Write(text4[0]) )
        self.play( Write(text4[1]) )
        self.wait(1)
        self.play( ReplacementTransform(text4,text4a) )
        self.wait(2)
        self.play( ReplacementTransform(text4a,text4b) )
        self.wait(4)
        self.play( ReplacementTransform(text4b,text4c) )
        self.wait(11)

        # determine dispersion relation
        self.play( Write(text5) )
        self.wait(2)
        self.play( ReplacementTransform(text5,text5a[0]) )
        self.wait(2)
        self.play( Write(text5a[1]) )
        self.wait(7)
        self.play( ReplacementTransform(text5a,text5b) )
        self.wait(1)
        self.play( ReplacementTransform(text5b,text5c) )
        self.wait(10)
        



class WaveText(Scene):
    def construct(self):
        title = Title(r"Wave speed is proportional to $\sqrt{\lambda}$")
        lst = write_lines(
            r"Can be deduced by",
            r"$\bullet$ experiment",
            r"$\bullet$ hydrodynamic model",
                r"incompressible, no viscosity,",
                r"does not feel bottom or edge",
            r"$\bullet$ dimensional analysis",
            indent = [0,0,0,1,1,0]
        ).shift(2*UP + 4*LEFT)
        framebox1 = SurroundingRectangle(lst[-1], buff = .1).set_color(BLUE)

        text2 = write_lines(
            r"Only remaining",
            r"parameters:",
            r"$\lambda \quad [ m ]$",
            r"$g \quad [ m/s^2 ]$",
            indent = [0,0,.5,.5]
        ).shift(.6*UP + 5*RIGHT).set_color(BLUE)
        brace2 = Brace( VGroup(lst[3],lst[4]), direction=RIGHT).set_color(BLUE)
        arrow2 = Arrow( brace2.get_corner(RIGHT), text2[1].get_corner(DOWN+LEFT) ).set_color(BLUE)

        text3 = write_lines(
            "Combination with units of velocity $[m/s]$:",
            r"$v = \sqrt{\lambda g}$",
            align_left=False
        ).shift(2.5*DOWN).set_color(BLUE)

        self.play(Write(title))
        self.wait(5)
        self.play(Write(lst[0]))
        self.wait(1)
        self.play(Write(lst[1]))
        self.wait(4)
        self.play(Write(lst[2]))
        self.wait(4)
        self.play(Write(lst[3]))
        self.wait(4)
        self.play(Write(lst[4]))
        self.wait(17)
        self.play(Write(lst[5]))
        self.wait()
        self.play(Create(framebox1))
        self.wait(5)
        self.play(Create(VGroup(brace2,arrow2)))
        self.play(Write(VGroup(*text2[0:2])))
        self.play(Write(text2[2]))
        self.wait(2)
        self.play(Write(text2[3]))
        self.wait(4)
        self.play(Write(text3[0]))
        self.wait(3)
        self.play(Write(text3[1]))
        self.wait(10)


class WavePacket(Scene):
    def construct(self):
        title = Title(r"Phase velocity and group velocity")
        ax1 = Axes(
            x_range=[-30, 30], y_range=[-2, 2], axis_config={"include_tip": False, "stroke_width": 0},
            y_length=8, x_length=16
        )
      
        t = ValueTracker(-2)
    
        ks = np.arange(.8,1.22,.01)
        kavg = 1.0
        func = lambda x,tt: sum([np.sin(k*x - omega(k)*tt) for k in ks])/len(ks)
        def get_plot1():
            lph = Line( 
                ax1.coords_to_point( vph(kavg)*t.get_value()-12 , 1 ),
                ax1.coords_to_point( vph(kavg)*t.get_value()-10 , 1 ),
                stroke_width=0
            )
            bph = Brace(lph, direction=[0,1,0], buff=0)
            bpht = bph.get_tex(r"v_{\mathrm{ph}}")
            textph = Tex("(phase velocity)").next_to(bpht, RIGHT, buff=.5)
            lgr = Line( 
                ax1.coords_to_point( vgr(kavg)*t.get_value()-8 , -1 ),
                ax1.coords_to_point( vgr(kavg)*t.get_value()+8 , -1 ),
                stroke_width=0
            )
            bgr = Brace(lgr, direction=[0,-1,0], buff=0)
            bgrt = bgr.get_tex(r"v_{\mathrm{gr}}")
            textgr = Tex("(group velocity)").next_to(bgrt, RIGHT, buff=.5)
            out = VGroup(
                ax1.plot(lambda x: func(x,t.get_value()), color=BLUE),
                lph,bph,bpht,textph,
                lgr,bgr,bgrt,textgr,
            )
            return out 
        graph1 = always_redraw( get_plot1 )

        plot1 = VGroup(ax1, graph1).shift(.5*DOWN)
        self.add(title)
        self.add(plot1)
        self.play( t.animate.set_value(5) , run_time=15, rate_func=linear )


class Superposition(Scene):
    def construct(self):
        title = Title(r"Superposition of two waves")

        # Axes for two sinewaves
        ax1 = Axes(
            x_range=[-24, 24], y_range=[-1.2, 1.2], axis_config={"include_tip": False},
            x_length=12, y_length=2,
        )
        labels1 = ax1.get_axis_labels(x_label="x", y_label=r"\sin(k x - \omega t)")
        ax2 = ax1.copy()
        # Reference frame for copies without axes shown
        ax1b = Axes(
            x_range=[-24, 24], y_range=[-1.2, 1.2], axis_config={"stroke_width": 0, "include_tip": False},
            x_length=12, y_length=2,
        )
        ax2b = ax1b.copy()
        labels2 = ax2.get_axis_labels(x_label="x", y_label=r"\sin(\bar k x - \bar \omega t)")
        # Axes for central superposition plot
        ax3 = Axes(
            x_range=[-24, 24], y_range=[-2.4, 2.4], axis_config={"include_tip": False},
            x_length=12, y_length=4,
        )
        labels3 = MathTex(r"\sin(k x - \omega t)+\sin(\bar k x - \bar \omega t)").shift(3*DOWN)

        t = ValueTracker(3)

        k1 = 1.0
        k2 = 1.2
        omega1 = omega(k1)
        omega2 = omega(k2)

        def get_plot1():
            return ax1.plot(lambda x: np.sin(k1*x - omega1*t.get_value()), color=BLUE)
        graph1 = always_redraw( get_plot1 )
        def get_plot2():
            return ax2.plot(lambda x: np.sin(k2*x - omega2*t.get_value()), color=BLUE)
        graph2 = always_redraw( get_plot2 )
        def get_plot1b():
            plot = ax1b.plot(lambda x: np.sin(k1*x - omega1*t.get_value()), color=BLUE, stroke_opacity=.5)
            return plot
        graph1b = always_redraw( get_plot1b )
        def get_plot2b():
            plot = ax2b.plot(lambda x: np.sin(k2*x - omega2*t.get_value()), color=BLUE, stroke_opacity=.5)
            return plot
        graph2b = always_redraw( get_plot2b )
        def get_plot3():
            return ax3.plot(lambda x: np.sin(k1*x - omega1*t.get_value()) + np.sin(k2*x - omega2*t.get_value()), color=WHITE)
        graph3 = always_redraw( get_plot3 )

        highlight3 = Circle().set_color(RED).shift(1*LEFT)
        highlight4 = Circle(radius=2.5).set_color(RED).shift(3*RIGHT)

        plot1 = VGroup(ax1, graph1, labels1).shift(1*UP+0*LEFT)
        plot2 = VGroup(ax2, graph2, labels2).shift(2*DOWN+0*LEFT)
        highlight1 = Circle(radius=.4).set_color(RED).move_to(labels1[1].get_center()).shift(.4*LEFT)
        highlight2 = Circle(radius=.4).set_color(RED).move_to(labels2[1].get_center()).shift(.4*LEFT)

        plot1b = VGroup(ax1b,graph1b).shift(UP+0*LEFT)
        plot2b = VGroup(ax2b,graph2b).shift(2*DOWN+0*LEFT)

        plot3 = VGroup(ax3, graph3, labels3)

        #Double graph
        self.add(title)
        self.wait()
        self.play(FadeIn(plot1, plot2))
        self.wait(5)
        self.play( ShowPassingFlash(highlight1, time_width=.5))
        self.play( ShowPassingFlash(highlight2, time_width=.5))
        #self.play( t.animate.set_value(3) , run_time=4 )
        self.wait()
        #self.play( t.animate.set_value(0) , run_time=.1 )

        #Transition
        self.play( FadeOut( ax1,ax2,labels1,labels2) )
        self.remove(graph1,graph2)
        self.add(plot1b,plot2b)
        #self.play( t.animate.set_value(3) , run_time=0 )
        self.play( plot1b.animate.move_to([0,0,0]), plot2b.animate.move_to([0,0,0]) )
        self.wait(2)

        #Single graph
        self.add(graph3)
        #self.play( t.animate.set_value(3) , run_time=0 )
        self.play( FadeIn(graph3,labels3) )
        self.wait()
        self.play( ShowPassingFlash(highlight4, time_width=.5))
        self.wait(3)
        self.play( ShowPassingFlash(highlight3, time_width=.5))
        self.wait(5)
        self.play( t.animate.set_value(15) , run_time=15, rate_func = linear)


class Superposition2(Scene):
    def construct(self):
        title = Title(r"Superposition of two waves") 
        # Axes for central superposition plot
        ax3 = Axes(
            x_range=[-24, 24], y_range=[-2.4, 2.4], axis_config={"include_tip": False},
            x_length=12, y_length=4,
        )

        t = ValueTracker(0)

        k1 = 1.0
        k2 = 1.2
        omega1 = omega(k1)
        omega2 = omega(k2)

        def get_plot3():
            return ax3.plot(lambda x: np.sin(k1*x - omega1*t.get_value()) + np.sin(k2*x - omega2*t.get_value()), color=WHITE)
        graph3 = always_redraw( get_plot3 )
        def get_plot4():
            plot = ax3.plot(lambda x: np.sin((k1+k2)/2*x - (omega1+omega2)/2*t.get_value()), color=BLUE) 
            return plot #DashedVMobject( plot, num_dashes=200)
        graph4 = always_redraw( get_plot4 )
        def get_plot5():
            plot = ax3.plot(lambda x: 2*np.cos((k1-k2)/2*x - (omega1-omega2)/2*t.get_value()) , color=YELLOW)
            return plot #DashedVMobject( plot , num_dashes=100 )
        graph5 = always_redraw( get_plot5 )

        text2 = VGroup(
            MathTex(r" \sin(k x - \omega t) + \sin(\bar k x - \bar \omega t)").move_to([0,2.2,0]),
            MathTex(r"=", r"\sin \left( \frac{k + \bar k}{2} x - \frac{\omega + \bar\omega}{2} t \right)",r" \cdot ",r"2 \cos \left( \frac{k - \bar k}{2} x - \frac{\omega - \bar\omega}{2}t \right)").move_to([0,-2.7,0]),
        )
        text2[0].set_color(WHITE)
        text2[1][1].set_color(BLUE)
        text2[1][3].set_color(YELLOW)

        self.add(title)
        self.add(graph3)
        self.wait()
        self.play(Write(text2[0]))
        self.wait(2)
        self.play(Create(graph4), run_time=2)
        self.wait(3)
        self.play(Create(graph5))
        self.play(Write(text2[1]))
        self.wait(4)
        self.play( t.animate.set_value(10) , run_time=16, rate_func = linear )




class SuperText(Scene):
    def construct(self):
        def write_lines2(content, align_left=True):
            out = [Tex(content[0])]
            halfline = 0
            indent = 1.7+out[-1].get_right()[0]
            for i in content[1:]:
                if type(i) == type(""):
                    if align_left:
                        out += [Tex(i).shift(.35*halfline*DOWN).align_to(out[0],LEFT).shift(indent*RIGHT)]
                        indent = 1.7+out[-1].get_right()[0]
                    else:
                        out += [Tex(i).shift(.35*halfline*DOWN)]
                if type(i) == type(1):
                    halfline += i
                    indent = 0
                if type(i) == type(.5):
                    indent = i
            return VGroup(*out)

        title = Title(r"Superposition of two waves")
        
        sinesum = MathTex(r" \sin",r"(k x - \omega t)",r" + \sin",r"(\bar k x - \bar \omega t)").set_color(WHITE).shift(2*UP+2*LEFT)
        sineprod = MathTex(r"= 2 ",r" \sin \left( \frac{k + \bar k}{2} x - \frac{\omega + \bar\omega}{2} t \right) ",r" \cos \left( \frac{k - \bar k}{2} x - \frac{\omega - \bar\omega}{2}t \right)").set_color(WHITE).align_to(sinesum,DOWN).shift(1.5*DOWN)

        trig = write_lines2(
            [r"Trig formulas:",
            r"$\sin(\alpha + \beta) = \sin\alpha \cos\beta + \cos\alpha\sin\beta$",
            2, 3.23, r"$\sin(\alpha - \beta) = \sin\alpha \cos\beta - \cos\alpha\sin\beta$",
            2, 2.0, r"$\Rightarrow$", r"$\sin(\alpha + \beta) + \sin(\alpha - \beta) = 2 \sin\alpha \cos\beta$",
            3, r"Use this identity with",
            2, 1.0, r"$\alpha = \frac{k + \bar k}{2} x - \frac{\omega + \bar\omega}{2}t$\quad and \quad$\beta  = \frac{k - \bar k}{2} x - \frac{\omega - \bar\omega}{2}t$",
            ]
        ).scale(.9).align_on_border(LEFT, buff=2).shift(.5*DOWN).set_color(GOLD)

        brace1 = Brace(sineprod[1], direction=[0,-1,0])
        brace2 = Brace(sineprod[2], direction=[0,-1,0])
        text1 = brace1.get_tex(r"\approx \sin(kx - \omega t)")
        text2 = write_lines2([
            r"wave with velocity", 
            3, r"$\displaystyle \frac{\frac{\omega - \bar\omega}{2} }{ \frac{k-\bar k}{2} }$"
        ],align_left=False).move_to(brace2.get_tip()).shift(1.2*DOWN)
        text2b = write_lines2([
            r"wave with velocity", 
            3, r"$\displaystyle\frac{\omega - \bar\omega}{ k-\bar k }$"
        ],align_left=False).move_to(brace2.get_tip()).shift(1.1*DOWN).set_color(YELLOW)
        text2b[1].shift(.1*UP)
        text2c = MathTex(r"\approx \frac{\mathrm{d} \omega}{\mathrm{d} k}").next_to(text2b[1],RIGHT).set_color(YELLOW)        
        text2d = MathTex(r"= v_{\mathrm{gr}}").next_to(text2c,RIGHT).shift(.1*DOWN).set_color(YELLOW)
        b1 = VGroup(brace1,text1).set_color(BLUE)
        b2 = VGroup(brace2,text2).set_color(YELLOW)
        
        highlight1 = Circle(radius=.8).set_color(RED).move_to(text2b[1].get_center())

        text3 = Tex(r"We found $\omega = \sqrt{360^\circ \cdot g \cdot k}$", r", so").align_on_border(LEFT, buff=.2).shift(2*DOWN)
        text4 = MathTex(r"v_{gr} = \frac{\mathrm{d} \omega}{\mathrm{d} k} = \frac12 \sqrt{\frac{360^\circ \cdot g}{k}} = \frac12 \frac{\omega}{k} = \frac12 v_{ph}").next_to(text3,DOWN,buff=.25).align_to(text3,LEFT)
        textvgr = VGroup(text3[0],text3[1],text4).set_color(GOLD).scale(.9)

        text5 = MathTex(r"v_{\mathrm{gr}}", r"= \frac12 \cdot", r"v_{\mathrm{ph}}").shift(2.5*DOWN + .8*LEFT)
        text6 = Tex(r"group velocity", r" = half of ", r"phase velocity").next_to(text5,DOWN,buff=.25)
        text5[2].set_color(BLUE)
        text5[0].set_color(YELLOW)
        text6[2].set_color(BLUE)
        text6[0].set_color(YELLOW)

        self.add(title)
        self.wait(3)
        self.play(Create(sinesum))
        self.wait(4)
        self.play(Create(sineprod), run_time=2)
        self.wait(1)
        self.play(FadeIn(trig))
        self.wait(4)
        self.play(FadeOut(trig))
        self.wait(4)
        self.play( Write(brace1) )
        self.wait(5)
        self.play( Write(text1) )
        self.wait(3)
        self.play( Write(brace2) )
        self.wait(11)
        self.play( Write(text2) )
        self.wait(2)
        self.play( ReplacementTransform(text2,text2b) )
        self.wait(45)
        self.play( ShowPassingFlash(highlight1,time_width=.5) )
        self.wait(2)
        self.play( Write(text2c) )
        self.wait(16)
        self.play( Write(text2d) )
        self.wait(2)
        self.play( FadeIn(textvgr[0]) )
        self.wait(6)
        self.play( FadeIn(textvgr[1:]) )
        self.wait(2)
        self.play(FadeOut(textvgr))
        self.wait(1)
        self.play(Create(text5))
        self.wait(2)
        self.play(Create(text6))
        self.wait(15)




