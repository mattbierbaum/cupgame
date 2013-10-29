The complex dynamics of beer pong
=================================

We all know that the game of beer pong is relatively simple in terms of
physics.  What this post presupposes is: maybe it isn't?

First, I feel like I have to explain what the game beer pong is.  And let me
tell you, it is exactly what it sounds like - a bunch of gentlemen getting
together in civil discourse, filling up solo cups with nice pure water, and
throwing ping pong balls into them in a rousy game of sport.  Well, fine, it's
dudes throwing balls into cups full of beer until their girlfriends leave for another party.

Now, what are typical questions that one may ask when thinking, "hmmm, I'm bored
and want to know about the physics of beer pong".  If you are asking yourself that 
question, there are likely other things in your life you should question. None the
less, you may arrive at these:

1. What is the precise trajectory I need to make it in a cup?
2. How accurate do I need to be to get the ball in the cup?  Something about angular error?
3. If I shot randomly at the cups what is the chance that I make it in? 

I have answers for these:

1. Easy, and useless.  You can't duplicate this with your hand. Give up. 
   Practice is better for this calculation.
2. Perform step 1, rinse and repeat.
3. A bit more subtle, but a good idea of its resolution comes from step 2, rinse, and sleep.

Let's ask a completely different question - what happens when we don't quite
make it?  What happens when the ball bounces off the rim?  Now there's
something to make your socks go up and down.

Okay, so everyone who has encountered physics knows that we can describe a ball
flying through air with a quadratic equation (except maybe Dwayne 'the rock'
Johnson who simply has no regard for gravity).  It looks something like this:

$$ y(t) = y_0 + v_{0,y}t - \frac{1}{2}g t^2 $$
$$ x(t) = x_0 + v_{0,x}t $$
$$ z(t) = z_0 + v_{0,z}t $$

where y is up, gravity is down and **there is no air resistance**.  I feel like
I have to bring up this point unncessarily early.  I do physics, not real
things.  Air resistance would make this real. Real hard.  So it's gone, like
Bruce Jenner's glory days. So, when you throw a ball, you can determine at any
time where it is, where it's going and which dreams it's left behind along the
way.  This is going to be useful to us later, so don't forget it.

Alright, so I guess we need to make some cups and balls and and then throw them
at one another.  Let's do it.

A cup rim is a torus.  A ball is a sphere. Bill Murray is my hero.  These are
just facts.  What is a little more contentious is that for the rest of this
post, I want to simplify the interaction of the cup rim and ball to a single
torus, leaving the ball to be an infinitessimal speck.  Imagine attaching one
of those awesome spherical neodymium magnets to a metal hoop and swinging it
around until the ball rolls over the entire surface.  The center of that sphere
would trace out the surface of a torus.  I hope that's convincing since I'm not
going to try again.  I tried to make a picture but it didn't turn out great.

** hilarious picture of this concept **

Finding collisions between this infinitessimal ball and ball-rim torus amounts
to is then finding the intersection of the equation of torus and a parabola. A
torus can be written

$$ ((x_n)^2 + (y_n)^2 + (z_n)^2 + R^2 - r^2)^2 = 4R^2(x_n^2 + y_n^2) $$

where $$x_n$$, $$y_n$$, $$z_n$$ are normalized coordinates, $$x_n = x - x_c$$
where $$x_c$$ is the center of the rim.  And all the rest too.  We can then
plug our previous equations and get an 8th order polynomial in time.  You can
convince yourself it's 8th order just by starring at the equation or realizing
how many times a cricket wicket could intersect a jelly doughnut on its side if
you were in a bitter lawsuit with Tim Horton's.

Okay, now that know every time (literally, the variable is t) the ball
intersects the rim, we have to figure out what happens when it bounces off.
It's just like bouncing on a table, except that table is very small and only
flat in a small region.  So, we find the normal of the rim at the collision
point, and reflect the ball about this normal.  Easy, breezy, misogynistic
cover girl.

Finally, we can start to look at some pictures a.k.a. did I just spew nonsense or 
is there a reason that this post goes on for 3 more pages. Let's drop ten
thousand ping pong balls from straight above a single cup rim and look at the
surface of interaction where they touch.  It looks like this:


The second picture is for an infinite array of cups, which we will be using throughout
Sections II-XVI.  I think this proves that I'm right so far.  So let's just go
for it.  I'll start again with ten thousands balls dropped directly overhead a
cup and color the pixel of where it started according to how many times the
ball bounced before it went into a cup.


