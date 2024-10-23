template <typename Jet>
Jet funkcja(const Jet & x, const Jet & y){
    Jet w =  sin( x*x - 2*(y+1))/exp(-y * y + cos(x * y));
    return w;
}
