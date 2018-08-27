import os
model_range=[]
model_range.append('WRF')
model_range.append('RAMS')

case_range=[]
case_range.append('CLN')
case_range.append('POL')

for model in model_range:
    for case in case_range:
        print('model: '+ model+ ' , case: '+ case)
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
                  +'"'
                  )
print('jobs submitted')

