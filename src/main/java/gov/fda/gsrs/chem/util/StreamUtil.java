package gov.fda.gsrs.chem.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.Supplier;
import java.util.stream.Collector;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import gov.nih.ncats.common.Tuple;
import gov.nih.ncats.common.util.CachedSupplier;

public class StreamUtil {
    

    public static <T> Stream<T> forIterator(Iterator<T> it){
        return gov.nih.ncats.common.stream.StreamUtil.forIterator(it);
    }

    /**
     * Creates a {@link Stream} for a given {@link Supplier},
     * similar to {@link Stream#generate(Supplier)}, except that
     * it is limited rather than infinite, and will return 
     * as soon as the provided {@link Supplier} returns an
     * empty optional.
     * @param sup
     * @return
     */
    public static <T> Stream<T> forGenerator(Supplier<Optional<T>> sup){
        Iterator<T> ir=new Iterator<T>(){
            public CachedSupplier<Optional<T>> next;
            public boolean initialized=false;

            @Override
            public synchronized boolean hasNext() {
                if(!initialized)initialize();
                return next.get().isPresent();
            }

            @Override
            public synchronized T next() {
                if(!initialized)initialize();
                Optional<T> n=next.get();
                cacheNext();
                return n.get();
            }

            public synchronized void cacheNext() {
                next.resetCache();
            }

            private void initialize(){
                next=CachedSupplier.of(sup);
                initialized=true;
            }
        };
        return forIterator(ir);
    }


    /**
     * Similar to {@link #forGenerator(Supplier)}, except that
     * it returning null triggers the termination of the stream,
     * rather than an empty {@link Optional}.
     * @param sup
     * @return
     */
    public static <T> Stream<T> forNullableGenerator(Supplier<T> sup){
        return forGenerator(()->{
            return Optional.ofNullable(sup.get());
        });
    }

    public static <T> Stream<T> forIterable(Iterable<T> sup){
        return StreamSupport.stream(sup.spliterator(),false);
    }

    private static <K,T> Stream<T> forNullableGenerator(final K k, Function<K,T> sup){
        return forNullableGenerator(()->sup.apply(k));
    }

    /**
     * Returns a {@link StreamGenerator}, which assumes that the provided
     * argument will be used to extract elements for a new stream. This is
     * useful for things which are like {@link Enumeration}s or {@link Iterator}
     * s, but do not explicitly implement those interfaces, making streams
     * slightly more difficult. A typical example may be something like this:
     * 
     * <pre>
     * <code>
     *     //Using standard state-full loop
     *     Matcher match = Pattern.compile("([A-Z])([0-9])").matcher("A1B2C3D4E5");
     *     StringBuilder sb = new StringBuilder();
     *     while(match.find()){
     *          if(sb.length()>0)sb.append(";");
     *          sb.append(match.group(1) + "." + match.group(2));
     *     }
     *     String mod1=sb.toString(); 
     *       
     *     
     *     //Using stream
     *     String mod2=StreamUtil.from(Pattern.compile("([A-Z])([0-9])").matcher("A1B2C3D4E5"))
     *          .streamWhile(m->m.find())
     *          .map(mr->mr.group(1) + "." + mr.group(2))
     *          .collect(Collectors.joining(";"));
     *      
     *     System.out.println(mod1);//A.1;B.2;C.3;D.4;E.5
     *     System.out.println(mod2);//A.1;B.2;C.3;D.4;E.5
     * </code>
     * </pre>
     * 
     * @param k
     * @return
     */

    public static <K> StreamGenerator<K> from(K k){
        
        
        return new StreamGenerator<K>(k);
    }


    public static class StreamGenerator<K>{
        private K k;
        private StreamGenerator(K k){
            this.k=k;
        }
        

        /**
         * Create a stream, where the seed value given is used
         * and function is provided to extract the "next" value
         * from that seed generator. If the function returns "null"
         * the stream will be terminated.
         * @param next
         * @return
         */
        public <T> Stream<T> streamNullable(ThrowableFunction<K,T> next){
            return forNullableGenerator(k, next);
        }

        /**
         * Create a stream, where the seed value given is used
         * and function is provided to extract the "next" value
         * from that seed generator. If the function returns an
         * empty Optional, then the stream is terminated.
         * @param next
         * @return
         */
        public <T> Stream<T> streamOptional(ThrowableFunction<K,Optional<T>> next){
            return forGenerator(()->next.apply(k));
        }

        /**
         * Create a stream, where the provided seed generator
         * is simply returned as long as a predicate test passes.
         * <p>
         * Note: the returned stream has undefined behavior
         * if parallelized prior to a mapping to a non-stateful
         * representation.
         * </p>
         * @param next
         * @return
         */
        public Stream<K> streamWhile(Predicate<K> next){
            return forGenerator(()->(next.test(k))?Optional.of(k):Optional.empty())
                    .sequential();
        }
    }

    /**
     * Utility function to help concatenate streams together
     * @param s
     * @return
     */
    public static <T> StreamConcatter<T> with(Stream<T> s){
        return new StreamConcatter<T>().and(s);
    }

    /**
     * Creates an infinite stream repeating the given
     * values in order.
     * @param elements
     * @return
     */
    public static <T> Stream<T> cycle(T ... elements){
        return Stream.generate(RolloverIterator.create(elements)::next);
    }


    private static final class RolloverIterator<T> {
        private final T[] elements;

        int index=0;
        private static <T> RolloverIterator<T> create(T[] elements){
            return new RolloverIterator<>(elements);
        }
        private  RolloverIterator(T[] elements){
            this.elements = elements;
        }


        public T next() {
            return elements[(index++)%elements.length];
        }
    }
    
    public static Stream<String> lines(String text){
        return Arrays.stream(text.split("\n"));
    }
    
    
    /**
     * Returns a supplier from the stream, which will return
     * the results of the stream, in order. Will return null
     * when the stream is exhausted.
     * @param stream
     * @return
     */
    public static <T> Supplier<T> supplierFor(Stream<T> stream){
        final Iterator<T> it = stream.iterator();
        Supplier<T> sup = ()->{
            synchronized(it){
                if(it.hasNext()){
                    return it.next();
                }
            }
            return null;
        };
        
        return sup;
    }


    /**
     * A simple Builder pattern for concatenating a
     * stream. This is just a convenience class, to
     * avoid having to call {@link Stream#concat(Stream, Stream)}
     * recursively on many items.
     * 
     * 
     * @author peryeata
     *
     * @param <T>
     */
    public static class StreamConcatter<T>{
        Stream<T> s= Stream.empty();
        private StreamConcatter(){}

        public StreamConcatter<T> and(Stream<T> newstream){
            s=Stream.concat(s, newstream);
            return this;
        }

        public StreamConcatter<T> and(Collection<T> newCollection){
            s=Stream.concat(s, newCollection.stream());
            return this;
        }

        public StreamConcatter<T> and(Iterable<T> newCollection){
            s=Stream.concat(s, forIterator(newCollection.iterator()));
            return this;
        }

        public StreamConcatter<T> and(T ... newThings){
            s=Stream.concat(s, Stream.of(newThings));
            return this;
        }
        
        public StreamConcatter<T> and(Supplier<T> newstream){
            s=Stream.concat(s, StreamUtil.forNullableGenerator(newstream));
            return this;
        }
        
        private StreamConcatter<T> andSt(Supplier<Stream<T>> newstream){
            s=Stream.concat(s, StreamUtil.forNullableGenerator(newstream).flatMap(o->o));
            return this;
        }
        /*
        public StreamConcatter<Tuple<Integer,T>> indexed(){
            return with(s.map(Util.toIndexedTuple()));
        }*/
        
        public Supplier<T> supplier(){
            return supplierFor(s);
        }
        
        public StreamConcatter<T> until(Predicate<T> pred){
            Supplier<T> sup=supplier();
            s=StreamUtil.forNullableGenerator(()->{
                T t= sup.get();
                if(t==null || pred.test(t)){
                    return null;
                }
                return t;
            });
            return this;
        }
        /**
         * chunk stream 
         * @param addIf
         * @return
         */
        public StreamConcatter<Stream<T>> chunk(long size){
           Final<CollectingStream<T>> cs = Final.of(new CollectingStream<T>(size));
           Stream<Stream<T>> meta=s.map(t->{
                                   CollectingStream<T> comb=cs.get();
                                   cs.set(comb.accept(t));
                                   if(cs.get()==comb)return Tuple.of(false,comb);
                                   return Tuple.of(true,comb);
                               })
                               .filter(t->t.k())
                               .map(c->c.v().stream());
           
           
           return with(meta).andSt(()->{
               return Optional.ofNullable(cs.getAndNullify())
                               .map(s->Stream.of(s.stream()))
                               .orElse(null);
           });
        }
        

        public Stream<T> stream(){
            return this.s;
        }
    }
    public static class Final<T>{
        private T t;
        public T get(){
            return t;
        }
        
        public T getAndNullify(){
            T t=get();
            set(null);
            return t;
        }
        public T getAndSet(T nt){
            T t=get();
            set(nt);
            return t;
        }
        public void set(T t){
            this.t=t;
        }
        public static<T>  Final<T> of(T t){
            Final<T> f=new Final<T>();
            f.set(t);
            return f;
        }
    }
    public static class CollectingStream<T>{
        List<T> list=new ArrayList<T>();
        final long max;
        
        public CollectingStream(long size) {
           max=size;
        }
        public synchronized CollectingStream<T> accept(T t){
            list.add(t);
            if(list.size()>=max){
                return new CollectingStream<T>(max);
            }
            return this;
        }
        public Stream<T> stream(){
            return list.stream();
        }
        
    }

    /**
     * This is just an extension of the {@link Function} interface
     * which will wrap a checked exception as a {@link RuntimeException}.
     * This is used for simplifying the Lambda functions that can throw an
     * error.
     * 
     */
    public static interface ThrowableFunction<T,V> extends Function<T,V>{

        @Override
        default V apply(T t){
            try {
                return applyThrowable(t);
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }

        public V applyThrowable(T t) throws Exception;
    }



    public static <T> Stream<T> forEnumeration(Enumeration<T> enumeration){
        return forNullableGenerator(()->{
            return (enumeration.hasMoreElements())?enumeration.nextElement():null;
        });
    }
    

    /**
     * <p>
     * Creates a {@link Collector} which will use the provided {@link Comparator}
     * and limit integer to collect only the "max" n elements, and return them
     * as a stream.
     * </p> 
     * 
     * <p>
     * For example, the following should give the same results:
     * </p>
     * 
     * <pre>
     *  <code>
     *      //Naive mechanism
     *      List<String> naive = Stream.of("B", "A", "E", "Z", "C", "Q", "T")
     *            .sort((a,b)->a.compareTo(b))
     *            .limit(3)
     *            .collect(Collectors.toList()); //["A", "B", "C"]
     *            
     *      //Targeted mechanism
     *      List<String> targeted = Stream.of("B", "A", "E", "Z", "C", "Q", "T")
     *            .collect(maxElements(3,(a,b)->a.compareTo(b))) //returns stream
     *            .collect(Collectors.toList()); //["A", "B", "C"]
     *   </code>
     * </pre>
     * 
     * @param n Limit of records returned
     * @param comp Comparator to do sorting
     * @return
     */
    public static <T> Collector<T,?,Stream<T>> maxElements(int n, Comparator<T> comp){
        return new ReducedCollector<T>(n,comp);
    }

    public static class ReducedCollector<T> implements Collector<T,TopNReducer<T>,Stream<T>>{
        private final Comparator<T> comp;
        private final int max;
        private ReducedCollector(int max,Comparator<T> comp){

            this.comp=comp;
            this.max=max;
        }

        @Override
        public Supplier<TopNReducer<T>> supplier() {
            return (()->new TopNReducer<T>(max,comp));
        }

        @Override
        public BiConsumer<TopNReducer<T>, T> accumulator() {

            return (red,t)->red.add(t);
        }

        @Override
        public BinaryOperator<TopNReducer<T>> combiner() {
            return (a,b)->{
                b.get().forEach(t->{
                    a.add(t);
                });
                return a;
            };
        }

        @Override
        public Function<TopNReducer<T>, Stream<T>> finisher() {
            return (red)->red.get();
        }

        @Override
        public Set<java.util.stream.Collector.Characteristics> characteristics() {
            return new HashSet<>();
        }


    }

    private static class TopNReducer<T> {
        private final PriorityQueue<T> pq;
        private final Comparator<T> comp;
        private final int cap;
        private final int effcap;
        private int _buff=0;


        public TopNReducer(int n, Comparator<T> comp) {
            this.comp=comp;
            pq = new PriorityQueue<T>(n, (a, b) -> {
                return -comp.compare(a, b);
            });
            cap = n;
            effcap = cap * 1;
        }

        public void add(T t) {
            pq.add(t);
            _buff++;
            if (_buff > effcap) {
                int r= _buff-cap;
                for(int i=0;i<r;i++){
                    pq.remove();
                }
                _buff=cap;
            }
        }

        public Stream<T> get() {
            return pq.stream()
                    .sorted(comp)
                    .limit(cap);
        }

    }
}