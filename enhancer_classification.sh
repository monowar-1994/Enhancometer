echo 'enhancer feature production Starting'
python3 tnc_with_pwm.py /home/rashik/Documents/Enhancometer/strong_enhancers.fasta /home/rashik/Documents/Enhancometer/enhancer_classification_data.csv 1
echo 'strong enhancer feature production complete'
python3 tnc_with_pwm.py /home/rashik/Documents/Enhancometer/weak_enhancers.fasta /home/rashik/Documents/Enhancometer/enhancer_classification_data.csv 1
echo 'weak enhancer feature production complete'
python3 tnc_with_pwm.py /home/rashik/Documents/Enhancometer/non_enhancers.fasta /home/rashik/Documents/Enhancometer/enhancer_classification_data.csv 0
echo 'non enhancer feature production complete'