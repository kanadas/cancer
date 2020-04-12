---
output: pdf_document
documentclass: extarticle
fontsize: 12pt
papersize: A4
geometry: "margin=1in"
header-includes: 
	- \pagenumbering{gobble}
---

# Optymalna strategia podawania leku

## Postawienie problemu

* układ równań, sformułowanie zadanie optymalizacyjnego
* dyskretyzacja układu równań, dyskretne zadanie optymalizacyjne

Prawdopodobnie dobrą strategią może być postawienie (a potem implementacja metod dla) zadania ciut ogólniejszego, niż to jedno konkretne?

## Cel minimum

* Dyskretyzacja za pomocą przybliżania sterowania (co najmniej) funkcją kawałkami stałą. Sprowadzenie problemu do zadania optymalizacji nieliniowej skończonego wymiaru.
* Implementacja rozwiązania w MATLAB-ie/Octave z wykorzystaniem gotowych narzędzi do optymalizacji nieliniowej (`Optimization Toolbox`).
	* wybór metody
	* uwzględnienie specyfiki tego konkretnego zadania (i MATLAB-a/Octave)
* Siatka jednorodna, modyfikowana do niejednorodnej przez użytkownika a posteriori.
* Znalezienie rozwiązania początkowego przy którym optymalizacja zbiega i daje możliwie dobry wynik.
* Testy numeryczne metody - weryfikacja (tego raczej nie będzie w pracy, ale trzeba to zrobic)
* Opis (o wyważonej objętości) zastosowanych metod i motywacja dokonanych wyborów.
* Rozwiązanie zadania dla kilku zestawów parametrów i dyskusja otrzymanych wyników
    * poprawność rozwiązań przy użyciu różnych metod
    * szybkość zbieżności/działania wykorzystanych metod
* Krytyka otrzymanych wyników
    * wskazanie potencjalnych problemów i konsekwencji z nich wynikających (klasy rozwiązań? nieciągłość funcji sterującej?)
    * dyskusja nt. ograniczeń użytych metod

## Rozszerzenia

* porównanie wyników z alternatywnymi metodami (biblioteki?)
* stabilność rozwiązań ze względu na parametry układu
* stabilność względem warunków początkowych
* Przybliżanie sterowania funkcją kawałkami wielomianową/splajnem.
* Implementacja algorytmu automatycznego zagęszczania siatki dyskretyzacji.
* Przybliżanie sterowania funkcją kawałkami wielomianową o zmiennym stopniu wielomianu i automatyczne znajdowanie odpowiedniego stopnia wielomianu na danym przedziale.
* Modyfikacja lub własna implementacja algorytmu optymalizacji nieliniowej uzwględniająca specyfikę zadania.
* Wykorzystanie specjalnych metod do walki z nieciągłością optymalnego sterowania.
* ...i co tylko jeszcze uda się zrobić ciekawego!
