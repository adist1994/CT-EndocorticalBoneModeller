/* 
 * File:   itkOptimisationCallback.h
 * Author: rap58
 *
 * Created on 14 November 2014, 12:11
 */

#ifndef ITKOPTIMISATIONCALLBACK_H
#define	ITKOPTIMISATIONCALLBACK_H

#include "itkCommand.h"

namespace itk
{
    template < typename TOptimizer >
    class OptimiserCallback : public Command {
    public:
        typedef OptimiserCallback             Self;
        typedef Command                  Superclass;
        typedef SmartPointer<Self>       Pointer;
        typedef SmartPointer<const Self> ConstPointer;
        itkTypeMacro( IterationCallback, Superclass );
        itkNewMacro( Self );
        typedef    TOptimizer     OptimizerType;
        void SetOptimizer( OptimizerType * optimizer ) {
            m_Optimizer = optimizer;
            observerTag = m_Optimizer->AddObserver( itk::IterationEvent(), this );
        }
        void TurnOffOptimizerObservation() {
            if(observerTag != -1) {
                m_Optimizer->RemoveObserver(observerTag);
                observerTag = -1;
            }
        }
        void TurnOnOptimizerObservation() {
            if ( m_Optimizer != ITK_NULLPTR && observerTag == -1) {
                observerTag = m_Optimizer->AddObserver( itk::IterationEvent(), this );
            }
        }
        void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE {
            Execute( (const itk::Object *)caller, event);
        }
        void Execute(const itk::Object *, const itk::EventObject & event) ITK_OVERRIDE {
            if( typeid( event ) == typeid( itk::StartEvent ) ) {
                std::clog << std::endl << "Position              Value" << std::endl << std::endl;
            } else if( typeid( event ) == typeid( itk::IterationEvent ) ) {
                //std::clog << m_Optimizer->GetCurrentIteration() << "   ";
                //std::clog << m_Optimizer->GetValue() << "   ";
                std::clog << m_Optimizer->GetCurrentPosition() << std::endl;
            } else if( typeid( event ) == typeid( itk::EndEvent ) ) {
                std::clog << std::endl << std::endl;
                //std::clog << "After " << m_Optimizer->GetCurrentIteration() << "  iterations " << std::endl;
                std::clog << "Solution is    = " << m_Optimizer->GetCurrentPosition() << std::endl;
            }
        }
        //  Software Guide : EndCodeSnippet
    protected:
        OptimiserCallback() {observerTag = -1; m_Optimizer = ITK_NULLPTR; };
        itk::WeakPointer<OptimizerType>   m_Optimizer;
        long observerTag;
    };
}
#endif	/* ITKOPTIMISATIONCALLBACK_H */

