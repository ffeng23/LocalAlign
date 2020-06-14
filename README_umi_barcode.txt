this is program used to extract umi and handling barcode (like cutadapt)

1)it is for umi_tools 
2)it takes umi pattern and extract, trim and add the umi to the sequence string. 
3)later we will expand this to cut adapters off the sequence like cutadapt.

2020-01-13

reason to do this because the umi_tools don't have anchor. To be precise, umi_tools can do "anchor", but the anchor is at the fixed location (nt positoin). So for example, let's say for some reason the barcode and sequence were shifted for one nt, then umi_tools can not find the anchor at the specific location and everything fell apart. But it is reasonable to allow the anchor to be "flexible" and use the achor to anchor the string and find the correct one if possible. 

(read more in help of ngs_umi_barcode -h)

see the folder ./UMI_barcode for testing and example.

