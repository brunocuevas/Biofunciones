package Biofunciones;

use strict;
use warnings;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION    =   1.00 ;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(&fasta_check &display_fasta_check);
@EXPORT_OK  =   (
                '&fasta_seq',  '&fasta_head',
                '&sequence_source',  '&type_seq',
                '&residues_proportion',  '&seq_len',
                '&display_residues_proportion',  '&dna_gc_calc',
                '&show_dna_gc', '&molecular_weigth_calculator',
                '&display_weigth', '&complementary_na',
                '&display_complementary', '&directorf',
                '&direct_indirectorf',
                '&pIcalculation', '&hidrophobicity'
                );
%EXPORT_TAGS=   ( );
our %aminoacid_hidrophobicity = (
    'A' =>1.800  ,
    'R' =>-4.500 ,
    'N' =>-3.500 ,
    'D' =>-3.500 ,
    'C' =>2.500  ,
    'Q' =>-3.500 ,
    'E' =>-3.500 ,
    'G' =>-0.400 ,
    'H' =>-3.200 ,
    'I' =>4.500  ,
    'L' =>3.800  ,
    'K' =>-3.900 ,
    'M' =>1.900  ,
    'F' =>2.800  ,
    'P' =>-1.600 ,
    'S' =>-0.800 ,
    'T' =>-0.700 ,
    'W' =>-0.900 ,
    'Y' =>-1.300 ,
    'V' =>4.200
); # Author(s): Kyte J., Doolittle R.F.
our %aminoacid_weights = (
    'A' =>  89,
    'R' =>  174,
    'N' =>  132,
    'D' =>  133,
    'C' =>  121,
    'Q' =>  146,
    'E' =>  147,
    'G' =>  75,
    'H' =>  155,
    'I' =>  131,
    'L' =>  131,
    'K' =>  146,
    'M' =>  149,
    'F' =>  165,
    'P' =>  115,
    'S' =>  105,
    'T' =>  119,
    'W' =>  204,
    'Y' =>  181,
    'V' =>  117
);
our %dna_weigths = (
    'A' => 331,
    'C' => 307,
    'G' => 342,
    'T' => 322
);
our %rna_weights = (
    'A' => 347,
    'C' => 323,
    'G' => 363,
    'U' => 324
);
our %codons = (
    'ATT'=>'I',
    'ATC'=>'I',
    'ATA'=>'I',
    'CTT'=>'L',
    'CTC'=>'L',
    'CTA'=>'L',
    'CTG'=>'L',
    'TTA'=>'L',
    'TTG'=>'L',
    'GTT'=>'V',
    'GTC'=>'V',
    'GTA'=>'V',
    'GTG'=>'V',
    'TTT'=>'F',
    'TTC'=>'F',
    'ATG'=>'M',
    'TGT'=>'C',
    'TGC'=>'C',
    'GCT'=>'A',
    'GCC'=>'A',
    'GCA'=>'A',
    'GCG'=>'A',
    'GGT'=>'G',
    'GGC'=>'G',
    'GGA'=>'G',
    'GGG'=>'G',
    'CCT'=>'P',
    'CCC'=>'P',
    'CCA'=>'P',
    'CCG'=>'P',
    'ACT'=>'T',
    'ACC'=>'T',
    'ACA'=>'T',
    'ACG'=>'T',
    'TCT'=>'S',
    'TCC'=>'S',
    'TCA'=>'S',
    'TCG'=>'S',
    'AGT'=>'S',
    'AGC'=>'S',
    'TAT'=>'Y',
    'TAC'=>'Y',
    'TGG'=>'W',
    'CAA'=>'Q',
    'CAG'=>'Q',
    'AAT'=>'N',
    'AAC'=>'N',
    'CAT'=>'H',
    'CAC'=>'H',
    'GAA'=>'E',
    'GAG'=>'E',
    'GAT'=>'D',
    'GAC'=>'D',
    'AAA'=>'K',
    'AAG'=>'K',
    'CGT'=>'R',
    'CGC'=>'R',
    'CGA'=>'R',
    'CGG'=>'R',
    'AGA'=>'R',
    'AGG'=>'R',
    'TAG'=>'-',
    'TGA'=>'-',
    'TAA'=>'-'
);
my %pkcarboxi   =   (
    'A'=> 3.1,
    'R'=> 3.1,
    'N'=> 3.1,
    'D'=> 3.1,
    'C'=> 3.1,
    'E'=> 3.1,
    'Q'=> 3.1,
    'G'=> 3.1,
    'H'=> 3.1,
    'I'=> 3.1,
    'L'=> 3.1,
    'K'=> 3.1,
    'M'=> 3.1,
    'F'=> 3.1,
    'P'=> 3.1,
    'S'=> 3.1,
    'T'=> 3.1,
    'W'=> 3.1,
    'Y'=> 3.1,
    'V'=> 3.1
);
my %pkamino     =   (
    'A'=> 8.0,
    'R'=> 8.0,
    'N'=> 8.0,
    'D'=> 8.0,
    'C'=> 8.0,
    'E'=> 8.0,
    'Q'=> 8.0,
    'G'=> 8.0,
    'H'=> 8.0,
    'I'=> 8.0,
    'L'=> 8.0,
    'K'=> 8.0,
    'M'=> 8.0,
    'F'=> 8.0,
    'P'=> 8.0,
    'S'=> 8.0,
    'T'=> 8.0,
    'W'=> 8.0,
    'Y'=> 8.0,
    'V'=> 8.0
);
my %pksidechainacid =   (
    'D' => 4.4,
    'C' => 8.5,
    'E' => 4.4,
    'Y' => 10.0
);
my %pksidechainbasic=   (
    'R' => 12,
    'H' => 6.5,
    'K' => 10.0
);

# our %pkcarboxi   =   (
#     'A'=> 2.35,
#     'R'=> 2.01,
#     'N'=> 2.02,
#     'D'=> 2.10,
#     'C'=> 2.05,
#     'E'=> 2.10,
#     'Q'=> 2.17,
#     'G'=> 2.35,
#     'H'=> 1.77,
#     'I'=> 2.32,
#     'L'=> 2.33,
#     'K'=> 2.18,
#     'M'=> 2.28,
#     'F'=> 2.58,
#     'P'=> 2.00,
#     'S'=> 2.21,
#     'T'=> 2.09,
#     'W'=> 2.38,
#     'Y'=> 2.20,
#     'V'=> 2.29
# );
# our %pkamino     =   (
#     'A'=> 9.87,
#     'R'=> 9.04,
#     'N'=> 8.80,
#     'D'=> 9.82,
#     'C'=> 10.25,
#     'E'=> 9.47,
#     'Q'=> 9.13,
#     'G'=> 9.78,
#     'H'=> 9.18,
#     'I'=> 9.76,
#     'L'=> 9.74,
#     'K'=> 8.95,
#     'M'=> 9.21,
#     'F'=> 9.24,
#     'P'=> 10.60,
#     'S'=> 9.15,
#     'T'=> 9.10,
#     'W'=> 9.39,
#     'Y'=> 9.11,
#     'V'=> 9.72
# );
# our %pksidechainacid =   (
#     'D' => 3.86,
#     'C' => 8.00,
#     'E' => 4.07,
#     'Y' => 10.07
# );
# our %pksidechainbasic=   (
#     'R' => 12.48,
#     'H' => 6.10,
#     'K' => 10.53
# );
sub fasta_check {
    # Comprueba que se le haya introducido un fichero fasta al vector.
    #   Devuelve "0" si no se trata de un fasta, y un "1" si se trata de un
    #   fasta.

    my $exitstatus=1;
    my @fasta;

    if (@_) {

        @fasta = @_ ;

        if ($fasta[0] =~ /^\>.*\n/i) {

            shift(@fasta);

            for (@fasta) {

                ($_ =~ /^\>.*\n/i) and $exitstatus = 0;

            }

            if ($exitstatus == 1) {

                return $exitstatus;

            } else {

                return $exitstatus;
            }
        }
    }
    else {
        die "No information to analyse";
    }
}
sub display_fasta_check {
    if (@_){
        (&fasta_check(@_)) and print "\tFASTA file recognised";
        (&fasta_check(@_)) or print "\tNot FASTA file recognised";
    } else {
        die 'No arguments';
    }
}

sub fasta_seq {
    if (@_){
        if (&fasta_check(@_)) {

            my @fastafile = @_;
            shift @fastafile;
            for (@fastafile){
                chomp $_;
            }
            my $line=join('',@fastafile);
            return $line;
        } else {
            die "No fasta format"
        }
    } else {
        die "No information to analyse";
    }
}

sub fasta_head {
    if (@_){
        if (&fasta_check(@_)) {

            my $head = $_[0] ;
            $head=~ s/\n//i ;
            $head=~ s/^\>//i ;
            return $head ;

        } else{

            die "No fasta format";

        }
    } else {
        die "No information to analyse";
    }
}
sub sequence_source {
    #Detecta el tipo de secuencia con la que estamos trabajando. Determina el
    #   siguiente índice de valores:
    #       Secuencias no biologicas    0
    #       Secuencias indeterminadas   1
    #       Secuencias de ADN           2
    #       Secuencias de ARN           3
    #       Secuencias de proteinas     4
    if (@_){

        if (&fasta_check(@_)==1){

            my $seq=fasta_seq(@_);
            my $letter;
            my %sequence_pieces_hash = (
                "A"=>1,"B"=>0, "C"=>1, "D" => 4, "E"=> 4, "F"=>4, "G"=>1, "H"=>4,
                "I"=>4, "J" => 0, "K"=>4, "L"=>4, "M"=>4, "N"=>4, "P"=>4, "Q"=>4, "R"=>4,
                "S"=>4, "T"=>2, "U"=>3, "V"=>4,"W"=>4,"X"=>0, "Y"=>4, "Z"=>0
            );
            my $punctuation = 0;
            for (my $iter=0; $iter<length($seq); $iter ++){
                $letter = substr($seq, $iter, 1);
                ($sequence_pieces_hash{$letter}==0) and return 0;
                if ($punctuation < $sequence_pieces_hash{$letter}){
                    $punctuation = $sequence_pieces_hash{$letter} ;
                }
            }
            if ($punctuation==3){
                do{
                    $letter = substr ($seq,0,1);
                    ($sequence_pieces_hash{$letter} == 2) and return "1";
                    $seq =~ s/^.//;
                } while ($seq);
            }
            return $punctuation;
        }
    }
}
sub type_seq {
    if (@_) {
        if (&fasta_check(@_) == 1){

            my $indicator = &sequence_source(@_);
            ($indicator==0) and
                print "\n\tNon-biological-sense sequence";
            ($indicator==1) and
                print "\n\tUndeterminated Biological-sense sequence";
            ($indicator==2) and
                print "\n\tDNA sequence";
            ($indicator==3) and
                print "\n\tRNA sequence";
            ($indicator==4) and
                print "\n\tProtein sequence";

        }
    }
}
sub residues_proportion {
    # Establece un hash con los valores de frecuencia de cada residuo
    #   dentro de la cadena.
    if (@_) {

        if (&fasta_check(@_)){

            my $seq = &fasta_seq(@_);
            my %residues_hash;
            my $letter;

            for (my $iter=0; $iter<length($seq) ; $iter ++ ){

                $letter = substr($seq, $iter, 1);
                if (exists $residues_hash{$letter}) {

                    $residues_hash{$letter} ++;

                } else {

                    $residues_hash{$letter} = 1;

                }
            }
            return %residues_hash;
        } else {
            die 'No fasta format';
        }
    }
}
sub seq_len {
    # Proporciona un valor de longitud de la secuencia.
    if (@_){
        if (&fasta_check(@_)){

            my $seq = &fasta_seq(@_);
            return length($seq);

        } else {
            die 'No Fasta Format';
        }
    }
}

sub display_residues_proportion {
    # Muestra en pantalla los residuos que componen una secuencia, con
    #   sus valores de frecuencia absoluta y relativa.
    if (@_){
        if (&fasta_check(@_)){

            my %residues_hash = &residues_proportion(@_);
            my $iter = 0;
            my $len = &seq_len;
            my $freq;
            for (sort(keys %residues_hash)) {

                $freq = $residues_hash{$_} / $len;
                print "\t\t\[$iter\]\t$_\t$residues_hash{$_}\t$freq\n";
                $iter ++;

            }
            print "\t\tLENGTH\t$len\n";
        } else {
            die 'No Fasta Format'
        }
    }
}
sub dna_gc_calc {
    # Calcula los valores de G+C dentro de la cadena.
    if (@_) {
        if (&fasta_check(@_)) {

            if ((&sequence_source(@_)==2) or (&sequence_source(@_))==3){

                my %residues = &residues_proportion(@_);
                my $long = &seq_len(@_);
                my $gcvalue = $residues{"G"} + $residues{"C"};
                return ($gcvalue/$long);
            } else {
                die 'No NUCLEIC sequence';
            }
        } else {
            die 'No fasta seq';
        }
    }
}
sub show_dna_gc {
    # Muestra en pantalla el valor de G+C.
    if (@_){
        my $value = dna_gc_calc(@_);
        print "\tG+C\t=\t$value";
    }
}
sub molecular_weigth_calculator {
    if (@_) {
        if (&fasta_check(@_)){
            my @file=@_;
            my $typeseq = &sequence_source(@_);
            if (($typeseq != 0) and ($typeseq != 1)) {
                my $reverse_control=1;
                my $reverse_control_2=0;
                my @seq=@_;
                my $weigth = 0;
                while ($reverse_control){

                    $reverse_control_2==1 and $reverse_control=0;
                    my %composition_values=&residues_proportion(@file);

                    my %current_hash;
                    ($typeseq == 4 ) and %current_hash=%aminoacid_weights;
                    ($typeseq == 3 ) and %current_hash=%rna_weights;
                    ($typeseq == 2 ) and %current_hash=%dna_weigths;

                    for (keys %composition_values) {
                        $weigth += ($current_hash{$_}-18)*$composition_values{$_};
                    }
                    (($typeseq == 2) or ($typeseq == 3)) and $reverse_control_2=1;
                    ($reverse_control_2==1) and @seq=complementary_na(@file);
                    ($typeseq == 4) and $reverse_control=0;
                }

                return $weigth;
            } else {
				
                die 'No Biological-Meaning sequence';
            }
        } else {
            die 'Fasta format required';
        }
    } else {
        die 'No arguments';
    }
}
sub display_weigth {
    if (@_){
        print "\t\t", &molecular_weigth_calculator(@_);
        print "\t Da"
    } else {
        die 'No arguments';
    }
}
sub complementary_na {
    if (@_) {
        if (fasta_check(@_)) {
            if ((sequence_source(@_)==2) or (sequence_source(@_)==3)) {
                my $seq = fasta_seq(@_);
                $seq = reverse($seq);
                (sequence_source(@_)==2) and $seq =~ tr/AGCT/TCGA/;
                (sequence_source(@_)==3) and $seq =~ tr/AGCU/UGCA/;
                return $seq;
            } else {
                die 'NO Nucleic ACID';
            }
        } else {
            die 'No Fasta';
        }
    } else {
        die 'No arguments';
    }
}
sub display_complementary {
    if (@_){
        my $seq = &complementary_na(@_);
        print "\n\t$seq\n";
    } else {
        die 'No arguments';
    }
}
sub codon2aminoacid {
    if (@_){
        my $codon = $_[0];
        if (length($codon)==3) {
            for (keys %codons) {
                ($_ eq $codon) and return %codons{$_};
            }
        }else {
            die 'Not Codon';
        }
    } else {
        die 'No arguments';
    }
}
sub trimorf {
    if (@_){
        my $seq_to_trim     = $_[0];
        my $long            = length($seq_to_trim);
        my $placetotrim     = $long - int($long/3);
        for (my $iter       = 0 ; $iter<$placetotrim ; $iter ++) {
            chop $seq_to_trim;
        }
        return $seq_to_trim;
    }
}
sub direct_indirectorf {
    if (@_){
        if (&fasta_check(@_)){
            if (&sequence_source(@_)==2){
                my @ref_seq_array       =   @_  ;
                my @output                      ;
                my @indirect_ref_seq            ;
                my $nx                  =   ">\n";
                $indirect_ref_seq[0]       =   &complementary_na(@_);
                @indirect_ref_seq       =   ($nx, @indirect_ref_seq);
                @output                 =   &directorf(@ref_seq_array);
                @output                 =   (@output, &directorf(@indirect_ref_seq));
                return @output          ;
            }
        }
    }
}
sub directorf {
    if (@_) {
        if (fasta_check(@_)) {
            if (sequence_source(@_)==2) {

                my @ref_seq_array;
                my $seq_x1   =   fasta_seq(@_);
                my $seq_x2   =   $seq_x1;
                $seq_x2   =~  s/^.//i;
                my $seq_x3   =   $seq_x2;
                $seq_x3   =~  s/^.//i;

                @ref_seq_array   =   (
                        \$seq_x1,
                        \$seq_x2,
                        \$seq_x3
                );

                #my $codon_number;
                my @output_orf;
                my $iter_2 = 0;

                for (@ref_seq_array) {
                    $output_orf[$iter_2]=&translate($$_);
                    $iter_2 ++;


                }
                #print "\nTHis is OUTPUT ORF =\n@output_orf";
                return @output_orf;

            } else {
                die 'No DNA'
            }
        }
    } else {
        die 'No arguments'
    }
}

sub translate {
    if(@_){
        my $seq     =   $_[0];
        my $len     =   int(length($seq)/3);
        my $output  =   "";
        for (my $iter = 0; $iter<($len*3) ; $iter += 3) {
            $output =   $output.&codon2aminoacid(substr($seq, $iter, 3));
        }
        return $output;
    }
}

sub pIcalculation {
    if (@_){
        if (&fasta_check(@_) and &sequence_source==4) {
            my @fasta      =       @_;
            my $seq        =       &fasta_seq(@fasta);
            my @seq        =       split('', $seq );
            my @positives   ;
            my @negatives   ;
            my $amino       ;
            my %positives_sidechains;
            my %negatives_sidechains;
            my $charge  =   1;
            my $pH      =   0;
            my $iter    =   0;
            for (@seq) {

                $amino = $_;
                (exists $pksidechainbasic{$amino}) and unshift(@positives, $pksidechainbasic{$amino});
                (exists $pksidechainacid{$amino}) and unshift(@negatives, $pksidechainacid{$amino});


            }

            @positives      =   sort({ $a <=> $b } @positives);
            @negatives      =   sort({ $a <=> $b } @negatives);

            for (@positives) {
                ((exists $positives_sidechains{$_}) and ($positives_sidechains{$_} ++))
                    or ($positives_sidechains{$_} = 1);
            }
            for (@negatives) {
                ((exists $negatives_sidechains{$_}) and ($negatives_sidechains{$_} ++))
                    or ($negatives_sidechains{$_} = 1);
            }
            while ($charge >= 0) {
                $pH += 0.001 ;
                $charge = &qnet(\%positives_sidechains, \%negatives_sidechains, $pH);
                $charge += 1/(1+(10**($pH-$pkamino{$seq[0]})));
                $charge -= 1/(1+(10**($pkcarboxi{$seq[-1]}-$pH)));
            }
            return $pH;
        }

    }

}

sub qnet {
    ### Calculates charge at a given pH
    if (@_){
        my ($ref1, $ref2, $ref3)    =   @_;
        my %positive_hash           =   %$ref1;
        my %negative_hash           =   %$ref2;
        my $pH                      =   $ref3;
        my $positive                =   0;
        my $negative                =   0;

        for (keys %positive_hash) {
            $positive   +=   ($positive_hash{$_})/(1+(10**($pH-$_)));
        }
        for (keys %negative_hash) {
            $negative   +=   ($negative_hash{$_})/(1+(10**($_-$pH)));
        }

        return $positive-$negative;

    }
}

sub hidrophobicity {
    if (@_) {
        if (sequence_source(@_) == 4 ) {
			my @seq				= @_;
            my %composition         ;
            	%composition        = &residues_proportion(@_)	;
			my $hidrophobicity_score = 0						;
			for (keys %composition) {
				$hidrophobicity_score += $aminoacid_hidrophobicity{$_}*$composition{$_};
			}
			$hidrophobicity_score = $hidrophobicity_score/&seq_len(@seq);
			return $hidrophobicity_score;
        } else {
			die 'No prot seq';
		}
    } else {
		die 'No arguments';
	}
}

1;
