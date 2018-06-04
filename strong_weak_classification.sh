echo 'Strong vs weak feature production Starting'
python3 tnc_with_pwm.py /home/rashik/Documents/Enhancometer/strong_enhancers.fasta /home/rashik/Documents/Enhancometer/strong_weak_classification_data.txt 1
echo 'strong enhancer feature production complete'
python3 tnc_with_pwm.py /home/rashik/Documents/Enhancometer/weak_enhancers.fasta /home/rashik/Documents/Enhancometer/strong_weak_classification_data.txt 2
echo 'weak enhancer feature production complete'
