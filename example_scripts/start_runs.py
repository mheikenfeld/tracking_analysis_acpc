import os
#v_max_range=[3,5,10,20,30]
w_max_threshold_range=[[1,3,5,10]]
v_max_range=[30]
#w_max_threshold_range=[5]

model_range=[]
model_range.append('WRF')
model_range.append('RAMS')

case_range=[]
case_range.append('CLN')
case_range.append('POL')

for model in model_range:
    for case in case_range:
        for v_max in v_max_range:
            for w_max_threshold in w_max_threshold_range:
                print('model: '+ model+ ' , case: '+ case + ' , v_max: '+ str(v_max) +' , w_max_threshold: '+ str(w_max_threshold))
                os.system('bsub' +' '
                          +'-q short-serial' + ' '
                          +'-n 1'+ ' '
                          +'-R "rusage[mem=50000]" -M 50000'+ ' '
                          +'-o bsub_%J.out'+' '
                          +'-e bsub_%J.err'+' '
                          +'-W 12:00'+' '
                          +'-B' + ' '
                          +'"'
                          + 'sh run_Tracking.sh'+ ' '                           
                          + model + ' ' 
                          + case + ' ' 
                          + str(v_max) + ' ' 
                          + str(w_max_threshold)+ ' ' 
                          +'"'
                          )
print('jobs submitted')

