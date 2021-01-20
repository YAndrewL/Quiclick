# -*- coding: utf-8 -*-
# @Time    : 2020/8/24 9:16 上午
# @Author  : Yufan A. Liu

# 文件配套和打包：quiclick的执行文件，sprint执行文件，repeat.txt文件包，utilities工具包,bwa工具包,samtools工具包
# import rpy2
import os
import sys
import warnings
import multiprocessing
import time
import pandas as pd
from pathlib import Path


# 参数文档
def helpDoc():
    print("")
    print(
        "----------------------------------------------------------------------------------------------------------------")
    print(
        "                                          Users help documentation                                              ")
    print(
        "----------------------------------------------------------------------------------------------------------------")
    print("")
    print(
        " Parameters         Description                                 Default & Note                    Optional(Y/N) ")
    print(
        "    -w           Pathway of raw data                             Current path                           Y       ")
    print(
        "    -i           Pathway of index files                 Create automatically (ENSEMBL version)          Y       ")
    print(
        "    -m        Multiprocessing for mapping                Single thread (supported files < 6)            Y       ")
    print(
        "    -p          Threads for mapping                       16 (CANNOT use with -m altogether)            Y       ")
    print(
        "    -E          RNA editing analysis model               Carry out RNA editing analysis                 Y       ")
    print(
        "    -D  Whether to execute differential expression analysis             NO                              Y        ")
    print(
        "    -k      Whether to save genome index files                     Do not save                          Y       ")
    print(
        "    -s                Species                               N/A (necessary parameter)                   N       ")

    print("")
    print(
        "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print(
        "NOTICE: only human, mouse, zebrafish, celegans are supported                                                    ")
    print(
        "Contact liuyf19@mails.tsinghua.edu.cn when questions arise.                                                     ")
    print(
        "----------------------------------------------------------------------------------------------------------------")
    print("")
    sys.exit(0)


# Download references from ENSEMBL 在没有参考基因组的情况下进行创建
def indexBuilder(path, spe):
    def refsDownload(fastaPath, gtfPath):
        if os.system('ls *.fa.gz') == 512:  # 如果参考基因组不存在，下载参考基因组
            os.system('wget ' + fastaPath)
        elif os.system('gzip -t $(ls *.fa.gz)') != 0:
            os.system('wget ' + fastaPath)
        else:
            pass

        if os.system('ls *.gtf.gz'):  # 如果gtf文件不存在，下载gtf文件
            os.system('wget ' + gtfPath)
        elif os.system('gzip -t $(ls *.gtf.gz)') != 0:
            os.system('wget ' + gtfPath)
        else:
            pass
        return

    os.chdir(path)
    refsFile = Path("./refs")
    if not refsFile.is_dir():
        os.system('mkdir refs')
    os.chdir('./refs')  # 进入refs目录

    if spe == "human":
        faPath = "ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        annoPath = "ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz"
        refsDownload(faPath, annoPath)
    elif spe == "mouse":
        faPath = "ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
        annoPath = "ftp://ftp.ensembl.org/pub/release-101/gtf/mus_musculus/Mus_musculus.GRCm38.101.gtf.gz"
        refsDownload(faPath, annoPath)
    elif spe == "zebrafish":
        faPath = "ftp://ftp.ensembl.org/pub/release-101/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"
        annoPath = "ftp://ftp.ensembl.org/pub/release-101/gtf/danio_rerio/Danio_rerio.GRCz11.101.gtf.gz"
        refsDownload(faPath, annoPath)
    elif spe == "celegans":
        faPath = "ftp://ftp.ensembl.org/pub/release-101/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz"
        annoPath = "ftp://ftp.ensembl.org/pub/release-101/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.101.gtf.gz"
        refsDownload(faPath, annoPath)
    else:
        warnings.warn("Species value is not supported! ")
        helpDoc()

    if os.system('ls *.fa.gz') == 512 or os.system('gzip -t $(ls *.fa.gz)') != 0:
        raise Exception("Download FASTA Files Error! Please download manually.")
    elif os.system('ls *.gtf.gz') == 512 or os.system('gzip -t $(ls *.gtf.gz)') != 0:
        raise Exception("Download GTF Files Error! Please download manually.")
    # 解压参考基因组和gtf文件
    os.system('gunzip $(ls *.fa.gz)')
    os.system('gunzip $(ls *.gtf.gz)')
    fasta = "".join(os.popen('ls *.fa').read().split())
    gtf = "".join(os.popen('ls *.gtf').read().split())
    os.environ['fasta'] = str(fasta)
    os.environ['gtf'] = str(gtf)

    os.system(
        executeFilePath + 'STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ./ '
                          '--genomeFastaFiles $fasta --sjdbGTFfile $gtf --sjdbOverhang 139')  # 建立genome index
    return


# Trim process & function
def trimProcess(sampleName):
    os.chdir(sampleName)
    if os.system('ls *1.fq.gz') == 512 or os.system('ls *2.fq.gz') == 512:
        raise Exception("At least ONE FASTQ File does not exist.")
    else:
        fq1 = "".join(os.popen('ls *1.fq.gz').read().split())
        fq2 = "".join(os.popen('ls *2.fq.gz').read().split())
        os.environ['fq1'] = str(fq1)
        os.environ['fq2'] = str(fq2)
        os.system(
            executeFilePath + 'trim_galore --quality 25 --phred33 --stringency 3 --length 20 --paired -trim1 --gzip '
                              '--output_dir ./trimmed_data $fq1 $fq2 2> trimProcess.log')
    return


# Mapping process & function
def mappingProcess(SampleName, indexPath, threads):
    os.chdir(indexPath)
    fasta = "".join(os.popen('ls *.fa').read().split())
    gtf = "".join(os.popen('ls *.gtf').read().split())
    os.chdir(workPath)
    os.chdir('./' + SampleName + '/trimmed_data')
    if os.system('ls *1_val_1.fq.gz') == 512 or os.system('ls *2_val_2.fq.gz') == 512:
        warnings.warn("NO valid clean data, trim process will be executed.")
        trimProcess(SampleName)
    else:
        fq1 = "".join(os.popen('ls *1_val_1.fq.gz').read().split())
        fq2 = "".join(os.popen('ls *2_val_2.fq.gz').read().split())
        os.environ['fq1'] = str(fq1)
        os.environ['fq2'] = str(fq2)
        os.chdir('../')
        os.system('mkdir mapping_output')
        try:
            os.system(
                executeFilePath + 'STAR --runThreadN ' + str(
                    threads) + ' --genomeDir ' + indexPath +
                ' --sjdbOverhang 139 --readFilesIn'
                ' ./trimmed_data/$fq1 ./trimmed_data/$fq2 --sjdbGTFfile ' + indexPath + gtf +
                ' --readFilesCommand zcat --outFileNamePrefix'
                ' ./mapping_output/ --outSAMtype BAM SortedByCoordinate 2> mappingProcess.log'
            )
        except ValueError:
            print("Mapping process ERROR! Please try again.")
    return


def readsCountProcess(SampleName, indexPath):
    os.chdir(indexPath)
    gtf = "".join(os.popen('ls *.gtf').read().split())
    os.environ['gtf'] = str(gtf)
    os.chdir(workPath)
    os.chdir('./' + SampleName)
    os.system(executeFilePath + 'htseq-count -s no -f bam ./mapping_output/Aligned.sortedByCoord.out.bam '
              + indexPath + gtf + ' > ' + SampleName + '.count 2> readsCountProcess.log')


def editingIndexBuilder(path):
    os.chdir(path)
    fasta = "".join(os.popen('ls *.fa').read().split())
    gtf = "".join(os.popen('ls *.gtf').read().split())
    global fastaStore
    fastaStore = fasta
    print("Build RNA editing index files......")
    os.chdir(workPath)
    command = executeFilePath + 'sprint pp -t ' + path + gtf + ' ' + path + fasta + ' ' + executeFilePath + 'bwa 2> editingIndex.log'
    os.system(command)


def editingProcess(SampleName, indexPath, spe):
    os.chdir(indexPath)
    os.chdir(workPath)
    os.chdir(SampleName)
    os.chdir("./trimmed_data")
    fq1 = "".join(os.popen('ls *1_val_1.fq.gz').read().split())
    fq2 = "".join(os.popen('ls *2_val_2.fq.gz').read().split())
    os.system('gunzip -c ' + fq1 + ' > ' + SampleName + '_1_val_1.fq')
    os.system('gunzip -c ' + fq2 + ' > ' + SampleName + '_2_val_2.fq')
    fq1 = "".join(os.popen('ls *1_val_1.fq').read().split())
    fq2 = "".join(os.popen('ls *2_val_2.fq').read().split())
    # execute sprint
    os.chdir(workPath)
    os.chdir(SampleName)
    os.chdir("./trimmed_data")
    os.system(
        executeFilePath + 'sprint ma -rp ' + executeFilePath + 'repeatFiles/' + spe + '.txt' + ' -c 6 -p 12 -1 ' + fq1 + ' -2 ' + fq2 + ' ' +
        indexPath + fastaStore + ' ../RES_res ' + executeFilePath + 'bwa ' + executeFilePath + 'samtools 2> editingProcess.log')


def sprintAnnotator(SampleName, indexPath):
    print("sprintAnnotator is running for editing annotation.")
    os.chdir(indexPath)
    gtf = "".join(os.popen('ls *.gtf').read().split())
    # gtfFile整理
    gtfFile = pd.read_table(gtf, sep="\t", skiprows=5, header=None, low_memory=False)
    gtfFile = gtfFile[gtfFile.iloc[:, 2].str.contains('gene')]
    gtfRaw = gtfFile.iloc[:, [0, 3, 4, 6]]
    gtfRaw.columns = ['ChromType', 'Start', 'End', 'StrandType']
    gtfSplit = pd.DataFrame(gtfFile.iloc[:, [8]].iloc[:, 0].str.split(r";|\"", expand=True)).iloc[:, [1, 7, 10, 13]]
    gtfSplit.columns = ['GeneID', 'GeneName', 'DBSource', 'GeneBiotype']
    gtfBedFile = gtfRaw.join(gtfSplit).reset_index(drop=True)

    os.chdir(workPath)
    os.chdir(SampleName)
    os.chdir('./RES_res')
    # EditingTable读取和整理
    EditingTable = pd.read_table("SPRINT_identified_all.res", sep="\t", header=0, low_memory=False)
    EditingTable['AD'] = EditingTable['AD:DP'].str.split(":", expand=True)[0]
    EditingTable['DP'] = EditingTable['AD:DP'].str.split(":", expand=True)[1]
    EditingTable['Ratio'] = round(pd.to_numeric(EditingTable['AD']) / pd.to_numeric(EditingTable['DP']), 4)
    EditingTable = EditingTable.drop(['AD:DP', 'Start(0base)'], axis=1)
    EditingTable.columns = ['Chrom', 'BasePos', 'Type', 'SupportingReads', 'Strand', 'AD', 'DP', 'Ratio']

    getNameList = ['unknown'] * len(EditingTable)
    EditingTable['GeneName'] = getNameList
    geneNames = list(gtfBedFile['GeneName'])
    geneEditCount = pd.Series(0, index=geneNames)
    for gene in geneNames:
        genePos = geneNames.index(gene)
        chrom = str(gtfBedFile.ChromType[genePos])
        # chrom = 'chr' + str(gtfBedFile.ChromType[genePos])
        startPos = gtfBedFile.Start[genePos]
        endPos = gtfBedFile.End[genePos]
        sumNumber = len(EditingTable[(EditingTable['Chrom'] == chrom) & (EditingTable['BasePos'] >= startPos) & (
                    EditingTable['BasePos'] <= endPos)])
        geneEditCount[gene] = sumNumber
        print("\rAnnotation proceeding %s gene in total %s" % ((genePos + 1), len(geneNames)), end="")
        EditingTable.loc[
            (EditingTable.Chrom == chrom) & (EditingTable.BasePos > startPos) & (EditingTable.BasePos < endPos),
            'GeneName'] = gene
    print(
        "\nEditing annotator is already done! \nFile is saved as %s%s/RES_res/SPRINT_gene_annotation.res and %s%s/RES_res/Gene_edit_count.txt " % (
            workPath, SampleName, workPath, SampleName))
    EditingTable.to_csv("SPRINT_gene_annotation.res", sep='\t', index=0)
    geneEditCount.to_frame().to_csv("Gene_edit_count.txt", sep='\t', header=0)


def diffExpr(nameList):
    print("Please assign control and experimental group number, and separate by space.")
    for k in range(len(nameList)):
        print("%s     %s" % (k, nameList[k]))
    controlList = list(map(int, input("Control group number: ").split(" ")))
    experList = list(map(int, input("Experimental group number: ").split(" ")))
    if len(controlList) <= 1 or len(experList) <= 1:
        raise ValueError("Biological replicate is necessary! ")
    controlPathDict = {}
    experDict = {}
    for s in controlList:
        controlPathDict[nameList[s]] = workPath + nameList[s] + "/" + nameList[s] + ".count"
    for s in experList:
        experDict[nameList[s]] = workPath + nameList[s] + "/" + nameList[s] + ".count"

    rFile = open("Differential_Express_analysis.R", 'w')
    rFile.write(
        "library(\"DESeq2\")\n"
        "library(\"ggplot2\")\n"
    )
    for name in experDict:
        rFile.write(
            "EX" + str(list(experDict.keys()).index(name) + 1) + "<-read.table(\"" + experDict[
                name] + "\",head=FALSE,row.names=1,stringsAsFactors=F,check.names=F)\n"
        )
    for name in controlPathDict:
        rFile.write(
            "CT" + str(list(controlPathDict.keys()).index(name) + 1) + "<-read.table(\"" + controlPathDict[
                name] + "\",head=FALSE,row.names=1,stringsAsFactors=F,check.names=F)\n"
        )
    varsDict = {}
    for j in experDict:
        varsDict[j] = ["EX" + str(list(experDict.keys()).index(j) + 1), "EX",
                       "Rep" + str(list(experDict.keys()).index(j) + 1),
                       "EX" + str(list(experDict.keys()).index(j) + 1)]
    for j in controlPathDict:
        varsDict[j] = ["CT" + str(list(controlPathDict.keys()).index(j) + 1), "CT",
                       "Rep" + str(list(controlPathDict.keys()).index(j) + 1),
                       "CT" + str(list(controlPathDict.keys()).index(j) + 1)]
    rFile.write(
        "count_raw_combine<-cbind("
    )
    for keys in varsDict:
        if list(varsDict.keys()).index(keys) + 1 == len(varsDict.keys()):
            rFile.write(varsDict[keys][3])
        else:
            rFile.write(varsDict[keys][3] + ",")
    rFile.write(")\n")
    rFile.write("colnames(count_raw_combine)<-c(")
    for keys in varsDict:
        if list(varsDict.keys()).index(keys) + 1 == len(varsDict.keys()):
            rFile.write("\"" + varsDict[keys][0] + "\"")
        else:
            rFile.write("\"" + varsDict[keys][0] + "\"" + ",")
    rFile.write(")\n")
    rFile.write(
        "count_combine<-count_raw_combine[!grepl(\"^_\", rownames(count_raw_combine)),]\n"
        "gene_id<-rownames(count_combine)\n"
        "samples<-colnames(count_combine)\n"
        "sample_detected_count<-apply(count_combine, 1, function(x)all(x >= 3))\n"
        "table(sample_detected_count)\n"
        "expressed_genes<-gene_id[sample_detected_count]\n"
        "write(expressed_genes, \"Expressed_genes.txt\")\n"
        "count_filter<-as.matrix(count_combine[expressed_genes,])\n"
    )
    rFile.write("colData<-data.frame(row.names = samples, GROUP = factor(c(")
    for keys in varsDict:
        if list(varsDict.keys()).index(keys) + 1 == len(varsDict.keys()):
            rFile.write("\"" + varsDict[keys][1] + "\"")
        else:
            rFile.write("\"" + varsDict[keys][1] + "\"" + ",")
    rFile.write(")), replicate = c(")
    for keys in varsDict:
        if list(varsDict.keys()).index(keys) + 1 == len(varsDict.keys()):
            rFile.write("\"" + varsDict[keys][2] + "\"")
        else:
            rFile.write("\"" + varsDict[keys][2] + "\"" + ",")
    rFile.write(
        "))\n"
        "dds<-DESeqDataSetFromMatrix(countData=count_filter, colData=colData, design=~GROUP)\n"
        "normData<-estimateSizeFactors(dds)\n"
        "normData_matrix<-counts(normData, normalized=TRUE)\n"
        "normData_matrix_log<-log2(normData_matrix)\n"
        "dds$GROUP<-relevel(dds$GROUP,ref=\"CT\")\n"
        "dds<-DESeq(dds)\n"
        "res<-results(dds)\n"
        "res_order<-res[order(res$log2FoldChange),]\n"
        "write.table(res_order, file=\"DESeq2_res_order.txt\",quote=FALSE,sep=\"\t\")\n"
    )
    rFile.write(
        "#volcano plot"
        "pdf(\"volcano_plot.pdf\")\n"
        "res<-read.csv(\"DESeq2_res_order.txt\",head=T,row.names=1,stringsAsFactors=F,check.names=F,sep=\"\t\")\n"
        "Tags<-ifelse(res$padj<0.05&abs(res$log2FoldChange)>1,\"both\","
        "ifelse(abs(res$log2FoldChange)>1,\"|LogFC|>1\",ifelse(res$padj<0.05,\"FDR<0.05\",\"Not Significant\")))\n"
        "res_new<-data.frame(res,Tags)\n"
        "ggplot(data=res_new,aes(x=log2FoldChange,y=-log10(padj),color=Tags,))+"
        "ggtitle(\"Volcano Plot\")+theme(plot.title = element_text(hjust=0.5))+"
        "geom_point(alpha=1,size =0.5)+xlim(-10,10)+"
        "scale_color_manual(name = \"\", values = c(\"#0096FF\", \"#999999\", \"#000000\", \"#9933FF\"),limits = c(\"both\", \"|LogFC|>1\", \"FDR<0.05\",\"Not Significant\"))\n"
        "dev.off()"
    )
    rFile.close()


if __name__ == '__main__':
    i = 0
    if len(sys.argv) < 2:
        helpDoc()

    # 参数声明
    executeFilePath = os.path.split(os.path.realpath(__file__))[0] + '/'
    workPath = os.getcwd()
    multiMapping = False
    mappingThread = 16
    rnaEditing = False
    indexFileSave = False
    diffExprAnalysis = False

    #  Parameter Setting
    while i < len(sys.argv):
        if sys.argv[i] == '-w':
            path = sys.argv[i + 1]
            os.chdir(path)
            workPath = path
        elif sys.argv[i] == '-i':
            genomeIndexPath = sys.argv[i + 1]
        elif sys.argv[i] == '-s':
            species = sys.argv[i + 1]
        elif sys.argv[i] == '-m':
            multiMapping = True
        elif sys.argv[i] == '-p':
            mappingThread = sys.argv[i + 1]
            if not mappingThread.isdigit():
                raise KeyError("Thread must be a integer! ")
        elif sys.argv[i] == '-E':
            rnaEditing = True
        elif sys.argv[i] == '-k':
            indexFileSave = True
        elif sys.argv[i] == '-D':
            diffExprAnalysis = True
        else:
            pass
        i += 1

    sampleList = os.popen('ls').read().split()
    '''
    if not diffExprAnalysis:
        diffExpr(sampleList)
    '''

    # 参数合法性判断
    if 'species' not in dir():
        warnings.warn("Parameter SPECIES (-s) is necessary.")
        helpDoc()
    else:
        if species == "human" or species == "mouse" or species == "zebrafish" or species == "celegans":
            pass
        else:
            warnings.warn("Species value is not supported!")
            helpDoc()

    if multiMapping and len(sampleList) > 6:
        warnings.warn("ONLY files < 6 will be supported.")
        helpDoc()
    elif '-m' in sys.argv and '-p' in sys.argv:
        raise ValueError("-m and -p CANNOT be used together! ")

    try:
        #  执行Trim Process
        if len(sampleList) > 0:
            print("Trim Process. This may take a long time......")
            for i in range(len(sampleList)):
                locals()['p' + str(i)] = multiprocessing.Process(target=trimProcess, args=(sampleList[i],))
                locals()['p' + str(i)].start()
                print("Running trim process on File %d......" % (i + 1))
            for i in range(len(sampleList)):
                locals()['p' + str(i)].join()
        else:
            raise Exception("No VALID Sample!")
        os.chdir(workPath)  # 每次执行结束之后回到开始的工作路径
        sum = 300  # 设置倒计时时间
        timeflush = 0.25  # 设置屏幕刷新的间隔时间
        for i in range(0, int(sum / timeflush)):
            TagList = ["\\", "|", "/", "—"]
            index = i % 4
            print("\rSettle and validate files..... {}".format(TagList[index]), end="")
            time.sleep(timeflush)
        print("\nTrim process finished! ")

        # 建立序列比对和RNA editing index
        if 'genomeIndexPath' not in dir():
            indexBuilder(workPath, species)
        os.chdir(workPath)

        #  执行Mapping Process

        if not multiMapping:
            if 'genomeIndexPath' in dir():
                for s in range(len(sampleList)):
                    print("Start mapping process on File %d" % (s + 1))
                    mappingProcess(sampleList[s], genomeIndexPath, mappingThread)
                    print("Mapping process %d finished." % (s + 1))
            else:
                for s in range(len(sampleList)):
                    print("Start mapping process on File %d" % (s + 1))
                    mappingProcess(sampleList[s], workPath + 'refs/', mappingThread)
                    print("Mapping process %d finished." % (s + 1))
            os.chdir(workPath)  # 每次执行结束之后回到开始的工作路径
        else:
            if 'genomeIndexPath' in dir():
                for s in range(len(sampleList)):
                    locals()['p' + str(s)] = multiprocessing.Process(target=mappingProcess,
                                                                     args=(
                                                                         sampleList[s], genomeIndexPath, mappingThread))
                    locals()['p' + str(s)].start()
                    print("Start mapping process on File %d" % (s + 1))
                for s in range(len(sampleList)):
                    locals()['p' + str(s)].join()
                print("Mapping process finished! ")
            else:
                for s in range(len(sampleList)):
                    locals()['p' + str(s)] = multiprocessing.Process(target=mappingProcess,
                                                                     args=(
                                                                         sampleList[s], workPath + 'refs/',
                                                                         mappingThread))
                    locals()['p' + str(s)].start()
                    print("Start mapping process on File %d" % (s + 1))
                for s in range(len(sampleList)):
                    locals()['p' + str(s)].join()
                print("Mapping process finished! ")
            os.chdir(workPath)  # 每次执行结束之后回到开始的工作路径

        #  执行Count Process
        if 'genomeIndexPath' in dir():
            for s in range(len(sampleList)):
                locals()['p' + str(s)] = multiprocessing.Process(target=readsCountProcess,
                                                                 args=(sampleList[s], genomeIndexPath))
                locals()['p' + str(s)].start()
                print("Running reads count process on File %d......" % (s + 1))
            for s in range(len(sampleList)):
                locals()['p' + str(s)].join()
            print("Reads count process finished! ")
        else:
            for s in range(len(sampleList)):
                locals()['p' + str(s)] = multiprocessing.Process(target=readsCountProcess,
                                                                 args=(sampleList[s], workPath +
                                                                       'refs/'))
                locals()['p' + str(s)].start()
                print("Running reads count process on File %d......" % (s + 1))
            for s in range(len(sampleList)):
                locals()['p' + str(s)].join()
            print("Reads count process finished! ")

        if not rnaEditing:
            pass
        else:
            os.chdir(workPath)
            if 'genomeIndexPath' not in dir():
                editingIndexBuilder(workPath + 'refs/')
            else:
                editingIndexBuilder(genomeIndexPath)
            # 执行editing process
            if 'genomeIndexPath' in dir():
                for s in range(len(sampleList)):
                    locals()['p' + str(s)] = multiprocessing.Process(target=editingProcess,
                                                                     args=(sampleList[s], genomeIndexPath, species))
                    locals()['p' + str(s)].start()
                    print("Running RNA editing process on File %d......" % (s + 1))
                for s in range(len(sampleList)):
                    locals()['p' + str(s)].join()
                print("RNA editing process finished! ")
            else:
                for s in range(len(sampleList)):
                    locals()['p' + str(s)] = multiprocessing.Process(target=editingProcess,
                                                                     args=(sampleList[s], workPath +
                                                                           'refs/', species))
                    locals()['p' + str(s)].start()
                    print("Running RNA editing process on File %d......" % (s + 1))
                for s in range(len(sampleList)):
                    locals()['p' + str(s)].join()
                print("RNA editing process finished! ")

        if not rnaEditing:
            pass
        else:
            for i in range(len(sampleList)):
                print("RNA editing annotation analysis on file %s " % (i + 1))
                if 'genomeIndexPath' not in dir():
                    sprintAnnotator(sampleList[i], workPath + 'refs/')
                else:
                    sprintAnnotator(sampleList[i], genomeIndexPath)

        if 'genomeIndexPath' not in dir():
            if not indexFileSave:
                os.chdir(workPath)
                os.system('rm -r ./refs')
        print("Congratulations! All processes done! ")
    except KeyboardInterrupt:

        os.chdir(workPath)
        if 'genomeIndexPath' not in dir():
            refsFile = Path("./refs")
            if refsFile.is_dir():
                os.system('rm -r ./refs')
        for i in range(len(sampleList)):
            os.chdir(sampleList[i])
            fq1 = "".join(os.popen('ls *1.fq.gz').read().split())
            fq2 = "".join(os.popen('ls *2.fq.gz').read().split())
            fileList = os.popen('ls').read().split()
            for j in range(len(fileList)):
                if fileList[j] != fq1 and fileList[j] != fq2:
                    os.system('rm -rf ' + fileList[j])
            os.chdir(workPath)
