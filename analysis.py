    
def run_gravity(steps = 10):
        
    """
    Runs gravity analysis.
    Note that the model should be built before
    calling this function.
    
    Keyword arguments:
    steps -- total number of analysis steps

    """
    
    ops.initialize()
    # Records the response of a number of nodes at every converged step
    ops.recorder('Node', '-file', 
                os.path.join('FGU_RC3DF_files','Gravity_Reactions.out'),
        '-time','-node', *list(range(1,5)), '-dof', *list(range(1,7)), 'reaction')

    # enforces the constraints using the transformation method
    ops.constraints('Transformation')

    # RCM numberer uses the reverse Cuthill-McKee scheme to order the matrix equations
    ops.numberer('RCM')

    # Construct a BandGeneralSOE linear system of equation
    ops.system('BandGeneral')

    # Uses the norm of the left hand side solution vector of the matrix equation to
    # determine if convergence has been reached
    ops.test('NormDispIncr', 1.0e-6, 100, 0, 2)

    # Uses the Newton-Raphson algorithm to solve the nonlinear residual equation
    ops.algorithm('Newton')

    # Uses LoadControl integrator object
    ops.integrator('LoadControl', 1/steps)

    # Constructs the Static Analysis object
    ops.analysis('Static')

    # Records the current state of the model
    ops.record()
    # Performs the analysis
    ops.analyze(steps)    
    
    print("Gravity analysis Done!")    

    


def run_modal(n_evs = 3):
    
    """
    Runs Modal analysis.
    Note that the model should be built before
    calling this function.

    Keyword arguments:
    n_evs -- number of eigenvalues

    """
    
    ops.initialize()

    # Records Eigenvector entries for Node 1,3 & 5 @ dof 1 
    ops.recorder('Node', '-file',
            os.path.join('FGU_RC3DF_files', 'ModalAnalysis_Node_EigenVectors_EigenVec1.out'),                 
                 '-node', *list(range(5,9)), '-dof', *[1,2], 'eigen 1')
    ops.recorder('Node', '-file',
            os.path.join('FGU_RC3DF_files', 'ModalAnalysis_Node_EigenVectors_EigenVec2.out'),
                 '-node', *list(range(5,9)), '-dof', *[1,2], 'eigen 2')

    ops.recorder('Node', '-file',
            os.path.join('FGU_RC3DF_files', 'ModalAnalysis_Node_EigenVectors_EigenVec3.out'),
                 '-node', *list(range(5,9)), '-dof', *[1,2], 'eigen 3')                 

    # Constructs a transformation constraint handler, 
    # which enforces the constraints using the transformation method.
    ops.constraints('Transformation')

    # Constructs a Plain degree-of-freedom numbering object
    # to provide the mapping between the degrees-of-freedom at
    # the nodes and the equation numbers.
    ops.numberer('Plain')

    # Construct a BandGeneralSOE linear system of equation object
    ops.system('BandGen')

    # Uses the norm of the left hand side solution vector of the matrix equation to
    # determine if convergence has been reached
    ops.test('NormDispIncr', 1.0e-12, 25, 0, 2)

    # Uses the Newton-Raphson algorithm to solve the nonlinear residual equation
    ops.algorithm('Newton')

    # Create a Newmark integrator.
    ops.integrator('Newmark', 0.5, 0.25)

    # Constructs the Transient Analysis object
    ops.analysis('Transient')

    # Eigenvalue analysis
    lamda = np.array(ops.eigen(n_evs))

    # Writing Eigenvalue information to file
    with open(os.path.join('FGU_RC3DF_files', 'ModalAnalysis_Node_EigenVectors_EigenVal.out'), "w") as eig_file:
        # Writing data to a file
        eig_file.write("lambda omega period frequency\n")
        for l in lamda:
            line_to_write = [l, l**0.5, 2*np.pi/(l**0.5), (l**0.5)/(2*np.pi)]
            eig_file.write('{:2.6e} {:2.6e} {:2.6e} {:2.6e}'.format(*line_to_write))
            eig_file.write('\n')


    # Record eigenvectors 
    ops.record()    
    
    print("Modal analysis Done!")
    
    

def run_pushover(m_1, steps = 10000, direction = 'X'):
    
    """
    Runs Pushover analysis.
    Note that the model should be built before
    calling this function. Also, Gravity analysis
    should be called afterwards.
    
    Keyword arguments:
    m_1 -- lumped mass value
    steps -- total number of analysis steps
    direction -- 'X' or 'Y' 
    """    
    phi = 1.0
    
    if direction == 'X':
        file_names = ['Pushover_Horizontal_ReactionsX.out',
                      'Pushover_Story_DisplacementX.out']
        loads = [[n,m_1*phi,0,0,0,0,0] for n in range(5, 9)]
        d_o_f = 1
    else:
        file_names = ['Pushover_Horizontal_ReactionsY.out',
              'Pushover_Story_DisplacementY.out']
        loads = [[n,0,m_1*phi,0,0,0,0] for n in range(5, 9)]
        d_o_f = 2
        
    ## Records the response of a number of nodes at every converged step
    # Global behaviour
    # records horizontal reactions of node 1 to 4
    ops.recorder('Node', '-file',
                 os.path.join('FGU_RC3DF_files',file_names[0]),
                 '-time','-node', *list(range(1,5)), '-dof', d_o_f, 'reaction')
    # records horizontal displacements of node 3 & 5
    ops.recorder('Node','-file',
                  os.path.join('FGU_RC3DF_files',file_names[1]),
                 '-time','-node', *list(range(5,9)), '-dof', d_o_f, 'disp')
    

    # Measure analysis duration
    tic = time.time()


    # Apply lateral load based on first mode shape in x direction (EC8-1)
    ops.pattern('Plain', 2, 1)
    [ops.load(*l) for l in loads];


    # Define step parameters
    step = +1.0e-05
    number_of_steps = steps

    # Constructs a transformation constraint handler, 
    # which enforces the constraints using the transformation method.
    ops.constraints('Transformation')

    # RCM numberer uses the reverse Cuthill-McKee scheme to order the matrix equations
    ops.numberer('RCM')

    # Construct a BandGeneralSOE linear system of equation object
    ops.system('BandGen')

    # Uses the norm of the left hand side solution vector of the matrix equation to
    # determine if convergence has been reached
    ops.test('NormDispIncr', 0.000001, 100)

    # Line search increases the effectiveness of the Newton method
    # when convergence is slow due to roughness of the residual.
    ops.algorithm('NewtonLineSearch',True, 0.8,
                  1000, 0.1, 10.0)

    # Displacement Control tries to determine the time step that will
    # result in a displacement increment for a particular degree-of-freedom
    # at a node to be a prescribed value.
    # Target node is 5 and dof is 1
    ops.integrator('DisplacementControl',5, d_o_f, step)

    # Constructs the Static Analysis object
    ops.analysis('Static')

    # Records the current state of the model
    ops.record()

    # Performs the analysis
    ops.analyze(number_of_steps)

    toc = time.time()
    ops.wipe()

    print('Pushover Analysis in {} Done in {:1.2f} seconds'.format(direction,(toc-tic)))


def run_time_history(direction = 'X', g_motion_id = 1, scaling_id = 1,
                     lamda = 1, acc_dir = 'FGU_RC3DF_files/acc_1.txt',alpha = 0.29):

    """
    Runs Time history analysis.
    Note that the model should be built before
    calling this function. Also, Gravity analysis
    should be called afterwards.
    
    Keyword arguments:
    direction -- 'X' or 'Y'
    g_motion_id -- Groundmotion id (in case you run many GMs, like in an IDA)
    scaling_id -- Scaling id (in case you run many GMs, like in an IDA)
    lamda -- Scaling of the GM
    acc_dir -- file directory of GM to read from
    alpha -- set to 1 for full time-history analysis 
    """
    
    if direction == 'X':
        file_names = ['TimeHistory_Horizontal_ReactionsX'+str(g_motion_id)+'.'+str(scaling_id)+'.out',
                      'TimeHistory_Story_DisplacementX'+str(g_motion_id)+'.'+str(scaling_id)+'.out',
                      'TimeHistory_Story_AccelerationX'+str(g_motion_id)+'.'+str(scaling_id)+'.out']
        omega = np.loadtxt(os.path.join('FGU_RC3DF_files', 'ModalAnalysis_Node_EigenVectors_EigenVal.out'),skiprows=1)[0,1]
        d_o_f = 1        
    else:
        file_names = ['TimeHistory_Horizontal_ReactionsY'+str(g_motion_id)+'.'+str(scaling_id)+'.out',
                      'TimeHistory_Story_DisplacementY'+str(g_motion_id)+'.'+str(scaling_id)+'.out',
                      'TimeHistory_Story_AccelerationY'+str(g_motion_id)+'.'+str(scaling_id)+'.out']
        omega = np.loadtxt(os.path.join('FGU_RC3DF_files', 'ModalAnalysis_Node_EigenVectors_EigenVal.out'),skiprows=1)[2,1]
        d_o_f = 2    
        

    ## Records the response of a number of nodes at every converged step
    # Global behaviour
    # records horizontal reactions of node 1 to 4
    ops.recorder('Node','-file',
                 os.path.join('FGU_RC3DF_files', file_names[0]),
                '-time','-node',*list(range(1,5)),'-dof',d_o_f,'reaction')
    # records horizontal displacements of node 5 to 8
    ops.recorder('Node','-file',
                os.path.join('FGU_RC3DF_files', file_names[1]),
                '-time','-node',*list(range(5,9)),'-dof',d_o_f,'disp')
    # records horizontal accelerations of node 5 to 8
    ops.recorder('Node','-file',
                os.path.join('FGU_RC3DF_files', file_names[2]),
                '-time','-node',*list(range(5,9)),'-dof',d_o_f,'accel')
 


    # Reading omega for obraining Rayleigh damping model
    xi = 0.05
    a_R, b_R = (0, 2*xi/omega)

    ## Analysis Parameters
    accelerogram = np.loadtxt(acc_dir)      # Loads accelerogram file
    dt = 0.02                               # Time-Step
    n_steps = len(accelerogram)             # Number of steps
    tol = 1.0e-6                            # prescribed tolerance
    max_iter = 5000                         # Maximum number of iterations per step




    # Uses the norm of the left hand side solution vector of the matrix equation to
    # determine if convergence has been reached
    ops.test('NormDispIncr', tol, max_iter,0,0)

    # RCM numberer uses the reverse Cuthill-McKee scheme to order the matrix equations
    ops.numberer('RCM')

    # Construct a BandGeneralSOE linear system of equation object
    ops.system('BandGen')

    # The relationship between load factor and time is input by the user as a 
    # series of discrete points
    ops.timeSeries('Path', 2, '-dt', dt, '-values', *accelerogram, '-factor', lamda)

    # allows the user to apply a uniform excitation to a model acting in a certain direction
    ops.pattern('UniformExcitation', 3, d_o_f,'-accel', 2)

    # Constructs a transformation constraint handler, 
    # which enforces the constraints using the transformation method.
    ops.constraints('Transformation')

    # Create a Newmark integrator.
    ops.integrator('Newmark', 0.5, 0.25)

    # assign damping to all previously-defined elements and nodes
    ops.rayleigh(a_R, b_R, 0.0, 0.0)

    # Introduces line search to the Newton algorithm to solve the nonlinear residual equation
    ops.algorithm('NewtonLineSearch',True,False,False,False,0.8,100,0.1,10.0)

    # Constructs the Transient Analysis object
    ops.analysis('Transient')

    # Measure analysis duration
    t = 0
    ok = 0
    print('Running Time-History analysis with lambda=',lamda)
    start_time = time.time()
    final_time = ops.getTime() + n_steps*dt
    dt_analysis = 0.1*dt

    while (ok == 0 and t <= alpha*final_time):

        ok = ops.analyze(1, dt_analysis)
        t = ops.getTime()    
        
    finish_time = time.time()
    
    if ok == 0:
        print('Time-History Analysis in {} Done in {:1.2f} seconds'.format(direction, finish_time-start_time))
    else:
        print('Time-History Analysis in {} Failed in {:1.2f} seconds'.format(direction, finish_time-start_time))
    
    ops.wipe()



def reset_analysis():
    """
    Resets the analysis by setting time to 0,
    removing the recorders and wiping the analysis.
    """    
    
    # Reset for next analysis case
    ##  Set the time in the Domain to zero
    ops.setTime(0.0)
    ## Set the loads constant in the domain
    ops.loadConst()
    ## Remove all recorder objects.
    ops.remove('recorders')
    ## destroy all components of the Analysis object
    ops.wipeAnalysis()
    
