<?xml version="1.0" encoding="UTF-8"?>
<!-- kfish2genexml.xsl -->
<!DOCTYPE xsl:stylesheet [
<!ENTITY nbsp   "&#160;">
]>

<xsl:stylesheet version="1.1" 
   xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
   xmlns:java="http://xml.apache.org/xalan/java"
   exclude-result-prefixes="java" 
   > 
<!-- xsl:include href="/templates/header.xsl" / --> 
<xsl:output method="html" doctype-public="-//W3C//DTD HTML 4.0 Transitional//EN"/>

<!-- some inline urls : others,dbxref via cgi -->
<xsl:param name="gmapurl">http://server7.eugenes.org:8091/gbrowse/cgi-bin/gbrowse/killifish2/?name=</xsl:param>
<!-- xsl:param name="gmapthumb"> </xsl:param -->
<xsl:param name="gmapthumb">http://server7.eugenes.org:8091/gbrowse/cgi-bin/gbrowse_img/killifish2/?w=300;t=Gene13r;q=</xsl:param>
<!-- thumb add? bestfish intron bwrna2 Gene12a -->

<!-- xsl:param name="goreport">http://eugenes.org/fbservlet/goreport/</xsl:param -->
<xsl:param name="goreport"> </xsl:param>
<xsl:param name="uniprot">http://www.uniprot.org/entry/</xsl:param>

<xsl:param name="ncbigni">http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&amp;cmd=Retrieve&amp;tool=euGenes&amp;list_uids=</xsl:param>
<xsl:param name="ncbiaa">http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&amp;db=Protein&amp;doptcmdl=GenPept&amp;tool=euGenes&amp;term=</xsl:param>

<xsl:param name="fish11url">/genepage/fish11xml/</xsl:param>

<xsl:param name="stkbkgene">http://www.ensembl.org/Gasterosteus_aculeatus/Gene/Summary?g=</xsl:param>
<xsl:param name="medakagene">http://www.ensembl.org/Oryzias_latipes/Gene/Summary?g=</xsl:param>
<xsl:param name="tetrgene">http://www.ensembl.org/Tetraodon_nigroviridis/Gene/Summary?g=</xsl:param>
<xsl:param name="tilapiagene">http://www.ensembl.org/Oreochromis_niloticus/Gene/Summary?g=</xsl:param>
<xsl:param name="zfishgene">http://www.ensembl.org/Danio_rerio/Gene/Summary?g=</xsl:param>
<xsl:param name="platygene">http://www.ensembl.org/Xiphophorus_maculatus/Gene/Summary?g=</xsl:param>
<xsl:param name="spotgargene">http://www.ensembl.org/Lepisosteus_oculatus/Gene/Summary?g=</xsl:param>
<xsl:param name="mayzebrgene">http://www.ensembl.org/Maylandia_zebra/Gene/Summary?g=</xsl:param>

<xsl:param name="db_xref">/cgi-bin/db_xref?</xsl:param>

<xsl:param name="search">/lucegene_arthropod/search</xsl:param>
<xsl:param name="searchthis">/lucegene_arthropod/search?q=kfish2genexml-all:</xsl:param>
<xsl:param name="thisurllk">/lucegene_arthropod/lookup?lib=kfish2genexml&amp;id=</xsl:param>
<xsl:param name="thisurl">/genepage/killifish/</xsl:param>
<xsl:param name="genesurl">/genepage/killifish/</xsl:param>
<!--  http://arthropods.eugenes.org/genepage/killifish/Funhe2EKm014548t1 -->

<xsl:param name="aaurl">/lucegene_arthropod/lookup?lib=kfish2seqs&amp;id=</xsl:param>
<xsl:param name="esurl">/lucegene_arthropod/lookup?lib=kfish2seqs&amp;id=</xsl:param>

<xsl:variable name="geneid">
<xsl:value-of select="//GeneSummary/BASIC_INFORMATION/GeneID"/>
</xsl:variable>

<xsl:template match="/">
	<html>
	<head>
	<title> <xsl:value-of select="//GeneSummary/Title"/> </title>
  <style>
  body { background-color: white; }
  A { color: blue;  font-weight: normal; }
  A:link { text-decoration: none;  }
  A:visited { text-decoration: none; color: green;  }
  A:hover { text-decoration: underline;  }
  
  .commonblock { 
    padding: 4pt 2pt 2pt 2pt;
    text-indent:  8pt; 
    font-family:  sans-serif; 
    font-size:    medium;
    font-weight:  bold;
    border-style: inset; border-width: thin; 
    color: black; background-color: #b8c8f8; 
    }
  
  table { 
    font-size: small; background-color: white;
    float: left;  width: 33%; vertical-align: top;
    }
     
  table.header { 
    font-size: small;  background-color: #f0f0ff; 
    border-style: inset; float: none;   width: 100%;  vertical-align: top;   
    }
    
  table.genepage {
    font-size: small; background-color: white; float: none; width: 100%;
    border-color: #201080; border-width: 2pt 0pt 1pt 0pt; border-style: double;
    }

  table.footer { 
    margin: 4pt 2pt 2pt 2pt;
    font-size: small;   
    float: left; width: 100%;
    border-color: #d03010; background-color: #d8e8f0;
  	border-width: 2pt 0pt 0pt 0pt; border-style: double;
    }
    
  table.xxxBASIC_INFORMATION { 
    font-size: small;  background-color: #f0f0ff; 
    border-style: inset; vertical-align: top;   
    }

  table.LOCATION {  font-size: small; width: 49%; float: left; }
  table.SIMILAR_GENES {  font-size: small; width: 49%;  }
  table.FUNCTION {  font-size: small;  width: 49%; }  
  table.xxSIMILAR_GENES { font-size: small; float: left; width: 99%; }
  
  </style>

  <script type="text/javascript" language="javascript">
  <!-- only works for one id/obj; want all ids handled   -->
    function toggleVis(obj)  {
        var el = document.getElementById(obj);
        if ( el.style.display != 'none' ) { el.style.display = 'none'; }
        else {  el.style.display = ''; }
    }
  </script>
  
  
	</head>

	<body id="body">
  <xsl:call-template name="wfleabase_header" />
<!-- put footer w/ in div rest so it floats together -->
<div id="GeneSummary" >
  <xsl:apply-templates/>
  <xsl:call-template name="wfleabase_footer" />
</div>
	</body>
	</html>
</xsl:template>

<xsl:template match="GeneSummary">
<table id="genepage" class="genepage" ><tr><td>
	<h1 align="left">
	<xsl:value-of select="Type"/> 
  for 
	<i>
	<xsl:value-of select="BASIC_INFORMATION/Species"/>
	</i> 
  <xsl:text>  </xsl:text>
  
  <xsl:choose>
  <xsl:when test='BASIC_INFORMATION/Symbol/text()'>
    <i><xsl:value-of select="BASIC_INFORMATION/Symbol"/></i>
  </xsl:when>
  <xsl:otherwise>
    <xsl:value-of select="BASIC_INFORMATION/GeneID"/>
  </xsl:otherwise>
  </xsl:choose>
		
	</h1>

  <xsl:apply-templates mode="first" select="BASIC_INFORMATION" />
  <xsl:apply-templates mode="first" select="EXPRESSION" />
  <xsl:apply-templates mode="first" select="GENE_PRODUCT" />
  <!-- <xsl:apply-templates mode="first" select="LOCATION" /> -->
  <!-- <xsl:apply-templates mode="first" select="SIMILAR_GENES" />  -->
  <!-- need to drop above from apply ... -->
  <xsl:apply-templates/>
  <br/>
</td></tr></table>
</xsl:template>


<xsl:template name="wfleabase_footer" > 
  <table id="footer"  class="footer" >
  <tr><td>
  <center>
  Developed at the <a href="http://iubio.bio.indiana.edu/gil/">Genome Informatics Lab</a>,
  Indiana University Biology Department.
  </center>
  </td></tr></table>
</xsl:template>

<xsl:template name="wfleabase_header" > 
<div id="header" >
  <table id="header"  class="header"  width="100%">
  <tr><td align="center" >
  <a href="/">euGenes
  <!-- <br/> <img src="icons/kfish2.gif" height="32" align="middle" /> -->
  </a>
  </td>
  
  <td align="center" >
     <a href="/blast/">BLAST</a>
  <!-- |  <a href="/biomart/">BioMart</a> -->
  <!-- |  <a href="/gbrowse/">GBrowse Maps</a> -->
  |  <a href="/docs/">Help</a>
  </td>
  
  <td align="left" valign="middle">
  <form action="{$search}">
  <input type="submit" value="Search:" /> 
  <input name="query" size="20" />
  <select name="lib">
  <option value="kfish2genexml">Gene pages</option>
  <option value="kfish2seqs">Sequences</option>
  <option value="webdocs">Web docs</option>
  </select>
  <a href="{$search}?help=1">[?]</a>
  </form>
  </td></tr></table>
</div>
</xsl:template>

<xsl:template name="commonhead" > 
  <tr class="commonblock"> 
  <td colspan="2" class="commonblock">
  <xsl:value-of select="name()"/>
  </td>
  </tr>
</xsl:template>

<xsl:template name="commonlabel" > 
  <td width="15%"><b><xsl:value-of select="name()"/> </b>:</td>
</xsl:template>

<xsl:template name="commonlabelname" > 
  <xsl:param name="cname"/>
  <td width="15%"><b><xsl:value-of select="$cname"/> </b>:</td>
</xsl:template>

<xsl:template  mode="commonfield"  match="GeneSummary/*/*"> 
  <tr valign="top"> 
  <xsl:call-template name="commonlabel" />
  <td> <xsl:apply-templates select="." /></td>
  </tr>
</xsl:template>

<xsl:template  mode="overfield"  match="GeneSummary/*/*">
  <tr><td>
  <b><xsl:value-of select="name()"/> </b> <br/>
  <xsl:apply-templates select="." />
  </td>
  </tr>
</xsl:template>


<!-- top level blocks -->
<xsl:template priority="1" name="toplevel" match="GeneSummary/*"> 
  <table id="{name()}" class="{name()}" width="100%" >
  <xsl:call-template name="commonhead" />
  <xsl:apply-templates mode="commonfield" />
	</table>
</xsl:template>

<xsl:template priority="6"  match="BASIC_INFORMATION"> </xsl:template>
<xsl:template priority="3"  mode="first"  match="BASIC_INFORMATION">
  <xsl:call-template name="toplevel" />
</xsl:template>
<!-- hide to save space.. -->
<xsl:template  mode="commonfield"  match="BASIC_INFORMATION/isoform">  </xsl:template>
<xsl:template  mode="commonfield"  match="BASIC_INFORMATION/Species">  </xsl:template>

<xsl:template priority="3" mode="first" match="OLD_BASIC_INFORMATION"> 
  <table id="{name()}" class="{name()}"  width="100%" >
  <xsl:call-template name="commonhead" />
  <tr><td colspan="2">
    <table cellpadding="4" border="0"><tr valign="top">
    <xsl:for-each select="*">
    <td><b><xsl:value-of select="name()"/> </b>: <xsl:apply-templates /> </td>
    </xsl:for-each>
</tr>
  </table>
  </td></tr> </table>
</xsl:template>

<xsl:template priority="6"  match="OLD_LOCATION"> </xsl:template>
<xsl:template mode="first"  match="OLD_LOCATION"> <xsl:call-template name="toplevel" /></xsl:template>
<xsl:template priority="5"  match="LOCATION"> <xsl:call-template name="toplevel" /></xsl:template>

<xsl:template priority="6"  match="EXPRESSION"> </xsl:template>
<xsl:template mode="first"  match="EXPRESSION"> <xsl:call-template name="toplevel" /></xsl:template>

<xsl:template priority="6"  match="GENE_PRODUCT"> </xsl:template>
<xsl:template mode="first"  match="GENE_PRODUCT"> <xsl:call-template name="toplevel" /></xsl:template>
<xsl:template mode="commonfield"  match="GENE_PRODUCT/NameSource">  </xsl:template>

<!-- 
<xsl:template priority="6" match="GENE_ONTOLOGY"> 
  <table id="{name()}" class="{name()}" width="100%" >
  <xsl:call-template name="commonhead" />
  <tr> <td> 
    <xsl:for-each select="*">
    <xsl:apply-templates /> 
    </xsl:for-each>
  </td></tr> </table>
</xsl:template>
 -->


<!-- field content handlers, links, etc. -->

<xsl:template match="Genome_map"> 
<xsl:variable name="val"><xsl:call-template name='clean_loc'/></xsl:variable>
 <a href="{$gmapurl}{$val}"><xsl:apply-templates/><br/> 
 <img border="1" src="{$gmapthumb}{$val}"/></a>
</xsl:template>

<xsl:template name="clean_loc"> 
  <xsl:choose>
  <xsl:when test='@url'><xsl:value-of select="@url"/></xsl:when>
  <xsl:when test='@id'><xsl:value-of select="@id"/></xsl:when>
  <!-- Scaffold0:2433342-2457267:+/Scaffold9868:2334079-2334121:.  all now end w/ :orient -->
  <xsl:when test="contains(text(),':+')"><xsl:value-of select="substring-before(text(),':+')"/></xsl:when>
  <xsl:when test="contains(text(),':-')"><xsl:value-of select="substring-before(text(),':-')"/></xsl:when>
  <xsl:when test="contains(text(),':.')"><xsl:value-of select="substring-before(text(),':.')"/></xsl:when>
  <xsl:when test="contains(text(),',')"><xsl:value-of select="substring-before(text(),',')"/></xsl:when>
  <xsl:when test="contains(text(),' ')"><xsl:value-of select="substring-before(text(),' ')"/></xsl:when>
  <!-- <xsl:when test="contains(text(),'/')"><xsl:value-of select="substring-before(text(),'/')"/></xsl:when> -->
  <xsl:otherwise><xsl:value-of select="text()"/></xsl:otherwise>
  </xsl:choose>
</xsl:template>

<!-- 
<xsl:template  match="Transcript"> 
<a href="/cgi-bin/transpage?{@id}"> <xsl:apply-templates/> </a>, 
</xsl:template>
 -->


<xsl:template match="exprval"> 
  <xsl:apply-templates/> <xsl:text> </xsl:text> <xsl:value-of select="@group"/>, 
  <!-- <xsl:apply-templates/> <xsl:text>, &nbsp;</xsl:text> -->
</xsl:template>

<xsl:template  match="Species"> <i><xsl:apply-templates/></i> </xsl:template>
<xsl:template  match="ref/Title"> <i><xsl:apply-templates/></i> </xsl:template>
<xsl:template  match="ref/Source"> <b><xsl:apply-templates/></b> </xsl:template>
<xsl:template  match="ref/Identifier"> <u><xsl:apply-templates/></u> </xsl:template>

<!-- xsl:template  match="goterm"> <a href="{$goreport}{@id}"> <xsl:apply-templates/> </a>, </xsl:template -->

<xsl:template match="url"> 
<xsl:variable name="val"><xsl:call-template name='clean_hval'/></xsl:variable>
<a href="{$val}"> <xsl:apply-templates/> </a>
</xsl:template>

<xsl:template match="db_xref"> 
<xsl:variable name="val"><xsl:call-template name='clean_hval'/></xsl:variable>
<xsl:variable name="dburl"> 
  <xsl:choose>
  <xsl:when test="starts-with($val,'Funhe')"><xsl:value-of select="$esurl"/></xsl:when>
  <xsl:when test="starts-with($val,'NP')"><xsl:value-of select="$ncbiaa"/></xsl:when>
  <xsl:otherwise><xsl:value-of select="$db_xref"/></xsl:otherwise>
  </xsl:choose>
</xsl:variable>
  <a href="{$dburl}{$val}"> <xsl:apply-templates/> </a>, 
</xsl:template>

<xsl:template name="clean_hval"> 
  <xsl:choose>
  <xsl:when test='@url'><xsl:value-of select="@url"/></xsl:when>
  <xsl:when test='@id'><xsl:value-of select="@id"/></xsl:when>
  <xsl:when test="contains(text(),',')"><xsl:value-of select="substring-before(text(),',')"/></xsl:when>
  <xsl:when test="contains(text(),' ')"><xsl:value-of select="substring-before(text(),' ')"/></xsl:when>
  <xsl:when test="contains(text(),'/')"><xsl:value-of select="substring-before(text(),'/')"/></xsl:when>
  <xsl:otherwise><xsl:value-of select="text()"/></xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template name="clean_hvallist">
  <xsl:choose>
  <xsl:when test="contains(text(),',')"><xsl:value-of select="translate(text(),',','+')"/></xsl:when>
  <xsl:otherwise><xsl:call-template name='clean_hval'/></xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="acc"> <!-- db_xref -->
<xsl:variable name="spp"><xsl:call-template name='clean_sppdb'/></xsl:variable>
<xsl:variable name="val"><xsl:call-template name='clean_sppid'/></xsl:variable>
<xsl:variable name="dburl"> 
  <xsl:choose>
<!-- // bad here w/ Funhe below .. spp=kfish2:val=Funhe
  <xsl:when test="starts-with($spp,'kfish2')"><xsl:value-of select="$kfishgene"/></xsl:when>
  <xsl:when test="starts-with($spp,'killifish')"><xsl:value-of select="$kfishgene"/></xsl:when>
 -->
  <xsl:when test="starts-with($spp,'medaka')"><xsl:value-of select="$medakagene"/></xsl:when>
  <xsl:when test="starts-with($spp,'stickleback')"><xsl:value-of select="$stkbkgene"/></xsl:when>
  <xsl:when test="starts-with($spp,'tetraodon')"><xsl:value-of select="$tetrgene"/></xsl:when>
  <xsl:when test="starts-with($spp,'tilapia')"><xsl:value-of select="$tilapiagene"/></xsl:when>
  <xsl:when test="starts-with($spp,'zfish')"><xsl:value-of select="$zfishgene"/></xsl:when>
  <xsl:when test="starts-with($spp,'zebrafish')"><xsl:value-of select="$zfishgene"/></xsl:when>
  <xsl:when test="starts-with($spp,'human')"><xsl:value-of select="$uniprot"/></xsl:when>
  <xsl:when test="starts-with($spp,'mayzebr')"><xsl:value-of select="$mayzebrgene"/></xsl:when>
  <xsl:when test="starts-with($spp,'platyfish')"><xsl:value-of select="$platygene"/></xsl:when>
  <xsl:when test="starts-with($spp,'spotgar')"><xsl:value-of select="$spotgargene"/></xsl:when>

  <xsl:when test="starts-with($val,'FISH11')"><xsl:value-of select="$fish11url"/></xsl:when>
  <xsl:when test="starts-with($val,'Funhe')"><xsl:value-of select="$thisurl"/></xsl:when>
  
  <xsl:otherwise></xsl:otherwise> 
  </xsl:choose>
</xsl:variable>
  <xsl:choose>
   <xsl:when test="$dburl = ''"> <xsl:apply-templates/> </xsl:when>
   <xsl:otherwise><a href="{$dburl}{$val}"> <xsl:apply-templates/> </a></xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template name="clean_sppid"> 
  <xsl:choose>
  <xsl:when test="contains(text(),':')"><xsl:value-of select="substring-after(text(),':')"/></xsl:when>
  <xsl:when test="contains(text(),'_NOT')"><xsl:value-of select="substring-after(text(),'_NOT')"/></xsl:when>
  <xsl:otherwise><xsl:value-of select="text()"/></xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template name="clean_sppdb"> 
  <xsl:choose>
  <xsl:when test="contains(text(),':')"><xsl:value-of select="substring-before(text(),':')"/></xsl:when>
  <xsl:when test="contains(text(),'_NOT')"><xsl:value-of select="substring-before(text(),'_NOT')"/></xsl:when>
  <xsl:otherwise><xsl:value-of select="text()"/></xsl:otherwise>
  </xsl:choose>
</xsl:template>


<xsl:template priority="4" match="SIMILAR_GENES"> 
  <table id="{name()}"  class="{name()}" >
  <xsl:call-template name="commonhead" />
  <tr><td> <xsl:apply-templates mode="simtable"/> </td></tr> 
	</table>
</xsl:template>

<xsl:template priority="6" match="OLD_SIMILAR_GENES"> </xsl:template>
<xsl:template priority="6" match="OLD_SIMILAR_GENES1"> 
  <table id="{name()}"  class="{name()}" >
  <xsl:call-template name="commonhead" />
  <tr><td>
  <table class="{name()}">  
		<tr>
		<xsl:for-each select="./Ortholog[1]/*">
		 <th align="left"><xsl:value-of select="name()"/></th>
		</xsl:for-each>
		</tr>
    <xsl:apply-templates mode="simtable"/>
  </table>
  </td></tr> 
	</table>
</xsl:template>


<!-- 
<xsl:template mode="simtable" match="Paralog1"> 
  <tr valign="top"> <td> <b>Paralog </b> 
  <xsl:for-each select="*"> <xsl:apply-templates mode="simtable" select="."/> </xsl:for-each>
  </td></tr> 
</xsl:template>

<xsl:template mode="simtable" match="Ortholog1"> 
  <tr valign="top"> <td> <b>Ortholog </b> 
  <xsl:for-each select="*"> <xsl:apply-templates mode="simtable" select="."/> </xsl:for-each>
  </td></tr> 
</xsl:template>
 -->

<xsl:template mode="simtable" match="Ortholog"> 
  <tr valign="top"> <td>  
  <xsl:for-each select="*"> <xsl:apply-templates mode="simtable" select="."/> </xsl:for-each>
  </td></tr> 
</xsl:template>

<xsl:template mode="simtable" match="align"> 
  &nbsp; <xsl:apply-templates select="."/>% aln/idn
</xsl:template>

<xsl:template mode="simtable" match="acc">  <xsl:apply-templates select="."/> </xsl:template>
<xsl:template mode="simtable" match="orclass"> 
 <b><xsl:apply-templates select="."/> </b> &nbsp; </xsl:template>


<xsl:template  mode="commonfield" priority="6" match="ADDITIONAL_INFORMATION/other"> 
  <xsl:apply-templates mode="otherfield" />
</xsl:template>

<xsl:template mode="list" match="db_xref"><xsl:apply-templates select="."/>, </xsl:template>

<xsl:template  mode="commonfield" priority="6" match="FUNCTION/Protein_domains">
  <tr><td>
  <b><xsl:value-of select="name()"/> </b> <br/>
  <xsl:apply-templates mode="list" select="." />
  </td></tr>
</xsl:template>


<xsl:template priority="10" mode="otherfield" match="Score"></xsl:template>

<xsl:template match="NCBI_genemodel/gnomon"> 
<xsl:variable name="val"><xsl:call-template name='clean_hval'/></xsl:variable>
<xsl:choose>
<xsl:when test="starts-with($val,'X')"><a href="{$db_xref}RefProt:{$val}"><xsl:apply-templates/></a></xsl:when>
<xsl:when test="starts-with($val,'N')"><a href="{$db_xref}RefProt:{$val}"><xsl:apply-templates/></a></xsl:when>
<xsl:otherwise><xsl:apply-templates/></xsl:otherwise>
</xsl:choose>
</xsl:template>


<xsl:template  mode="otherfield"  match="*">
  <xsl:variable name="name">
  <xsl:choose>
  <xsl:when test="@type"><xsl:value-of select="@type"/></xsl:when>
  <xsl:when test="name() = 'upkw'">Keywords</xsl:when>
  <xsl:when test="name() = 'upde'">Description</xsl:when>
  <xsl:when test="name() = 'upcc'">Comments</xsl:when>
  <xsl:when test="name() = 'nrbits'">Bit score</xsl:when>
  <xsl:when test="name() = 'ncbiov'">NCBI gene</xsl:when>
  <xsl:when test="name() = 'tilex'">Expression</xsl:when>
  <xsl:otherwise><xsl:value-of select="name()"/></xsl:otherwise>
  </xsl:choose>
  </xsl:variable>
  <tr> 
  <!-- <xsl:call-template name="commonlabel" /> -->
  <xsl:call-template name="commonlabelname" >
    <xsl:with-param name="cname" select="$name"/>
  </xsl:call-template>
  <td> <xsl:apply-templates select="." /></td>
  </tr>
</xsl:template>

<xsl:template  match="Paralog/acc"> <!--  Paralog/acc was paralog-->
  <a href="{$genesurl}{text()}"> <xsl:apply-templates/> </a> 
</xsl:template>

<xsl:template  match="FamilyDbxref"> 
<xsl:variable name="ogurl"><xsl:call-template name='clean_omclid'/></xsl:variable>
  <a href="{$fish11url}{$ogurl}"> <xsl:apply-templates/> </a> 
</xsl:template>
<xsl:template name="clean_omclid"> <!-- fixme at fish11xml: FISH11G_G26.s3 << subids needed -->
  <xsl:choose>
  <xsl:when test="contains(text(),'.')"><xsl:value-of select="substring-before(text(),'.')"/></xsl:when>
  <xsl:otherwise><xsl:value-of select="text()"/></xsl:otherwise>
  </xsl:choose>
</xsl:template>




<!-- 
<xsl:template  match="ncbi|ncbiov">
<xsl:variable name="val"><xsl:call-template name='clean_hvallist'/></xsl:variable>
  <a href="{$thisurl}{$val}"> <xsl:apply-templates/> </a>
</xsl:template>
 -->


<xsl:template priority="6"  match="SUMMARY"> 
  <table id="{name()}"  class="{name()}" width="100%" style="width:auto">
  <xsl:call-template name="commonhead" />
  <xsl:apply-templates mode="sumtable" />
	</table>
</xsl:template>

<xsl:template mode="sumtable" match="SUMMARY/*"> <!-- summary/text summary/html .. -->
  <tr><td> <xsl:apply-templates /></td></tr> 
</xsl:template>

<xsl:strip-space elements="protein"/>
<xsl:template  match="GeneSummary/GENE_PRODUCT/Sequence/protein">
<!-- this only works for one element id=name; need new id for each new entry -->
<a href="#" onclick="toggleVis('{name()}')">&gt; Protein</a>
size:<xsl:value-of select="string-length(translate(text(),' &#xA;',''))"/>
<div id="{name()}" style="display: none;">
<pre>&gt;<xsl:value-of select="$geneid"/>
<xsl:value-of select="translate(text(),' ','')"/>
</pre>
</div>
</xsl:template>

<!-- // no good need other js for tag id or what.
<xsl:template  match="GeneSummary/GENE_PRODUCT/Sequence/protein2">
<a href="#" onclick="toggleVis('{id()}')">&gt; Protein</a>
size:<xsl:value-of select="string-length(translate(text(),' &#xA;',''))"/>
<div id="{id()}" style="display: none;">
<pre>&gt;<xsl:value-of select="$geneid"/><br/>
<xsl:value-of select="translate(text(),' ','')"/>
</pre>
</div>
</xsl:template>
 -->

<!-- hide metadata pulled at top -->
<xsl:template priority="10" match="GeneSummary/Title"></xsl:template>
<xsl:template priority="10" match="GeneSummary/Source"></xsl:template>
<xsl:template priority="10" match="GeneSummary/Type"></xsl:template>

</xsl:stylesheet>
