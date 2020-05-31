#include <bits/stdc++.h>
#include <include/SR_Engine.cpp>
#include <include/vec.cpp>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/ioctl.h>
#include <poll.h>
#include <omp.h>
#define TESTS
using namespace std;
using namespace std::chrono;

constexpr bool tests{true};
constexpr bool Make_Replays{false};
constexpr int PIPE_READ{0},PIPE_WRITE{1};
constexpr int N{1};
constexpr bool Debug_AI{false},Timeout{false};
constexpr double FirstTurnTime{1*(Timeout?1:10)},TimeLimit{0.15*(Timeout?1:10)};
constexpr char Replay_Filename[]{"Replay.txt"};
default_random_engine generator(system_clock::now().time_since_epoch().count());

bool stop{false};//Global flag to stop all arena threads when SIGTERM is received

struct AI{
    int id,pid,outPipe,errPipe,inPipe;
    string name;
    inline void stop(){
        if(alive()){
            kill(pid,SIGTERM);
            int status;
            waitpid(pid,&status,0);//It is necessary to read the exit code for the process to stop
            if(!WIFEXITED(status)){//If not exited normally try to "kill -9" the process
                kill(pid,SIGKILL);
            }
        }
    }
    inline bool alive()const{
        return kill(pid,0)!=-1;//Check if process is still running
    }
    inline void Feed_Inputs(const string &inputs){
        if(write(inPipe,&inputs[0],inputs.size())!=inputs.size()){
            throw(5);
        }
    }
    inline ~AI(){
        close(errPipe);
        close(outPipe);
        close(inPipe);
        stop();
    }
};

void StartProcess(AI &Bot){
    int StdinPipe[2];
    int StdoutPipe[2];
    int StderrPipe[2];
    if(pipe(StdinPipe)<0){
        perror("allocating pipe for child input redirect");
    }
    if(pipe(StdoutPipe)<0){
        close(StdinPipe[PIPE_READ]);
        close(StdinPipe[PIPE_WRITE]);
        perror("allocating pipe for child output redirect");
    }
    if(pipe(StderrPipe)<0){
        close(StderrPipe[PIPE_READ]);
        close(StderrPipe[PIPE_WRITE]);
        perror("allocating pipe for child stderr redirect failed");
    }
    int nchild{fork()};
    if(nchild==0){//Child process
        if(dup2(StdinPipe[PIPE_READ],STDIN_FILENO)==-1){// redirect stdin
            perror("redirecting stdin");
            return;
        }
        if(dup2(StdoutPipe[PIPE_WRITE],STDOUT_FILENO)==-1){// redirect stdout
            perror("redirecting stdout");
            return;
        }
        if(dup2(StderrPipe[PIPE_WRITE],STDERR_FILENO)==-1){// redirect stderr
            perror("redirecting stderr");
            return;
        }
        close(StdinPipe[PIPE_READ]);
        close(StdinPipe[PIPE_WRITE]);
        close(StdoutPipe[PIPE_READ]);
        close(StdoutPipe[PIPE_WRITE]);
        close(StderrPipe[PIPE_READ]);
        close(StderrPipe[PIPE_WRITE]);
        execl(Bot.name.c_str(),Bot.name.c_str(),(char*)NULL);//(char*)Null is really important
        //If you get past the previous line its an error
        perror("exec of the child process");
    }
    else if(nchild>0){//Parent process
        close(StdinPipe[PIPE_READ]);//Parent does not read from stdin of child
        close(StdoutPipe[PIPE_WRITE]);//Parent does not write to stdout of child
        close(StderrPipe[PIPE_WRITE]);//Parent does not write to stderr of child
        Bot.inPipe=StdinPipe[PIPE_WRITE];
        Bot.outPipe=StdoutPipe[PIPE_READ];
        Bot.errPipe=StderrPipe[PIPE_READ];
        Bot.pid=nchild;
    }
    else{//failed to create child
        close(StdinPipe[PIPE_READ]);
        close(StdinPipe[PIPE_WRITE]);
        close(StdoutPipe[PIPE_READ]);
        close(StdoutPipe[PIPE_WRITE]);
        perror("Failed to create child process");
    }
}

inline string EmptyPipe(const int fd){
    int nbytes;
    if(ioctl(fd,FIONREAD,&nbytes)<0){
        throw(4);
    }
    string out;
    out.resize(nbytes);
    if(read(fd,&out[0],nbytes)<0){
        throw(4);
    }
    return out;
}

bool IsValidMove(const string &Move){
    stringstream ss(Move);
    string line;
    getline(ss,line);
    stringstream ss2(line);
    vec target;
    int thrust;
    if(!(ss2 >> target >> thrust)){
        return false;
    }
    return true;
}

string GetMove(AI &Bot,const int turn){
    pollfd outpoll{Bot.outPipe,POLLIN};
    time_point<system_clock> Start_Time{system_clock::now()};
    string out;
    while(static_cast<duration<double>>(system_clock::now()-Start_Time).count()<(turn==1?FirstTurnTime:TimeLimit)){
        double TimeLeft{(turn==1?FirstTurnTime:TimeLimit)-static_cast<duration<double>>(system_clock::now()-Start_Time).count()};
        if(poll(&outpoll,1,TimeLeft)){
            out+=EmptyPipe(Bot.outPipe);
            if(IsValidMove(out)){
                return out;
            }
        }
    }
    throw(1);
}

inline Map Generate_Map(default_random_engine &generator)noexcept{
    uniform_int_distribution<int> Map_Distrib(0,Maps.size()-1);
    Map C{Maps[Map_Distrib(generator)]};
    return C;
}

state Random_Initial_State(){
    state S;
    S.timeout=Timeout_Length;
    S.C=Generate_Map(generator);
    vec CP0_CP1{S.C[1]-S.C[0]};
    CP0_CP1.normalise();
    const float angle{static_cast<float>(Angle(S.C[1]-S.C[0]))};
    S.p=pod{S.C[0],vec{0,0},angle,1,0};
    return S;
}

action StringToAction(const string &mv_str,const pod &p){
    stringstream ss(mv_str);
    action mv;
    int thrust;
    vec target;
    ss >> target >> thrust;
    mv.thrust=thrust;
    double desired_angle{Angle(target-p.r)};
    const double angle_given{Closest_Angle(p.angle,desired_angle)};
    /*ASSERT(thrust>=0 && thrust<=200,"Invalid thrust in StringToAction()");
    ASSERT(abs(Angular_Distance(desired_angle,p.angle))<19,"Requested delta_angle is "+to_string(abs(Angular_Distance(desired_angle,p.angle)))+" in StringToAction()");*/
    /*if(abs(Angular_Distance(desired_angle,p.angle))>18+1e-3){
        cerr << "Warning: Delta angle : " << Angular_Distance(desired_angle,p.angle) << endl;
        cerr << "Desired: " << desired_angle << " Current: " << p.angle << " Given: " << angle_given << endl;
    }*/
    const double angle_dist{Angular_Distance(angle_given,p.angle)};
    mv.delta_angle=Angular_Distance(desired_angle,Angle_Sum(p.angle,angle_dist))<Angular_Distance(desired_angle,Angle_Sum(p.angle,-angle_dist))?angle_dist:-angle_dist;
    if(tests && Angular_Distance(desired_angle,Angle_Sum(p.angle,mv.delta_angle))>1){
        cerr << "Target to delta_angle conversion is dodgy " << setprecision(3) << p.angle  << " " << angle_given << " " << mv.delta_angle << endl;
    }
    if(tests && abs(mv.delta_angle)>18+1e-3){
        cerr << "Error: Delta angle : " << mv.delta_angle << endl;
        cerr << "Desired: " << desired_angle << " Current: " << p.angle << " Given: " << angle_given << endl;
    }
    return mv;
}

int Play_Game(array <string,N> &Bot_Names){
    state S{Random_Initial_State()};
    vector<AI> Bot(Bot_Names.size());
    for(int id=0;id<Bot_Names.size();++id){
        Bot[id].id=id;
        Bot[id].name=Bot_Names[id];
        StartProcess(Bot[id]);
    }
    int turn{0};
    //Feed first turn inputs
    for(AI &b:Bot){
        stringstream ss;
        ss << S.C.size()*3 << endl;
        for(int i=0;i<3;++i){
            for(int j=1;j<=S.C.size();++j){
                ss << S.C[j%S.C.size()] << endl;
            }
        }
        b.Feed_Inputs(ss.str());
        //cerr << ss.str();
    }
    while(++turn>0 && !stop){
        for(int id=0;id<Bot.size();++id){
            if(Bot[id].alive()){
                //Feed turn inputs
                stringstream ss;
                ss << (S.p.next-1+(S.p.lap+(S.p.next==0?1:0))*S.C.size()) << " " << round(S.p.r.x) << " " << round(S.p.r.y) << " " << round(S.p.v.x) << " " << round(S.p.v.y) << " " << round(S.p.angle) << endl;
                //cerr << ss.str();
                try{
                    Bot[id].Feed_Inputs(ss.str());
                    string Move{GetMove(Bot[id],turn)};
                    string err_str{EmptyPipe(Bot[id].errPipe)};
                    if(Debug_AI){
                        ofstream err_out("log.txt",ios::app);
                        err_out << err_str << endl;
                    }
                    //cerr << Move;
                    const action mv{StringToAction(Move,S.p)};
                    ofstream replay_os{Replay_Filename,ios::app};
                    if(Make_Replays && omp_get_thread_num()==0){
                        for(const vec &base:S.C){
                            replay_os << base << " ";
                        }
                        replay_os << endl;
                        replay_os << S.p.r << " " << S.p.angle << " " << S.p.v << " " << S.p.next << endl;
                        replay_os << mv << endl;
                    }
                    Simulate(S,mv);
                    if(S.game_over()){
                        return S.finished_race()?S.turn:-2;
                    }
                }
                catch(int ex){
                    if(ex==1){//Timeout
                        cerr << "Loss by Timeout of AI " << Bot[id].id << " name: " << Bot[id].name << endl;
                    }
                    else if(ex==3){
                        cerr << "Invalid move from AI " << Bot[id].id << " name: " << Bot[id].name << endl;
                    }
                    else if(ex==4){
                        cerr << "Error emptying pipe of AI " << Bot[id].name << endl;
                    }
                    else if(ex==5){
                        cerr << "AI " << Bot[id].name << " died before being able to give it inputs" << endl;
                    }
                    else{
                        cerr << "AI " << Bot[id].name << " stopped after throw int " << ex << endl;
                    }
                    return -1;
                }
            }
            else{
                return -2;
            }
        }
    }
    return -2;
}

void StopArena(const int signum){
    stop=true;
}

int main(int argc,char **argv){
    if(argc<2){
        cerr << "Program takes 1 input, the names of the AI to test" << endl;
        return 0;
    }
    int N_Threads{1};
    if(argc>=3){//Optional N_Threads parameter
        N_Threads=min(2*omp_get_num_procs(),max(1,atoi(argv[2])));
        cerr << "Running " << N_Threads << " arena threads" << endl;
    }
    array<string,N> Bot_Names;
    for(int i=0;i<N;++i){
        Bot_Names[i]=argv[i+1];
    }
    cout << "Testing AI " << Bot_Names[0] << endl;
    for(int i=0;i<Bot_Names.size();++i){
        ifstream Test{Bot_Names[i].c_str()};
        if(!Test){
            cerr << Bot_Names[i] << " couldn't be found" << endl;
            return 0;
        }
        Test.close();
    }
    signal(SIGTERM,StopArena);//Register SIGTERM signal handler so the arena can cleanup when you kill it
    signal(SIGPIPE,SIG_IGN);//Ignore SIGPIPE to avoid the arena crashing when an AI crashes
    int maps_done{0},maps_failed{0};
    double average_turns{0};
    #pragma omp parallel num_threads(N_Threads) shared(average_turns,maps_done,Bot_Names)
    while(!stop){
        int turns{Play_Game(Bot_Names)};
        if(turns>0){
            #pragma omp atomic
            average_turns+=(turns-average_turns)/(++maps_done);
        }
        else{
            ++maps_failed;
            //cerr << "Bot didn't complete map" << endl;
        }
        #pragma omp critical
        cerr << "Average turns: " << setprecision(5) << average_turns << " with " << maps_done << " maps done with " << setprecision(2) << 100.0*maps_failed/(maps_failed+maps_done) << "% maps failed" << endl;
    }
}